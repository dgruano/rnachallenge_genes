"""
workflow/scripts/download_tool_datasets.py
==========================================
Stage 0 — Download tool dataset FASTAs.

For each classifier tool listed in config/tool_sources.yaml, attempt to
download its published FASTA/sequence files using the configured strategy.

Strategies implemented:
  github       — GitHub Contents API (recursive, handles subdirs)
  direct_list  — Direct URL list; each URL saved with a derived filename
  sourceforge  — SourceForge JSON directory listing → direct downloads
  osf          — Open Science Framework storage API
  readme_dir   — Fetch README from genouest server, parse linked filenames,
                 then try to download each; falls back to common extensions
  web          — Single URL attempt; logs a manual-download instruction on failure

Failed downloads are logged as WARNING and skipped — the pipeline continues.
All results are recorded in the output manifest JSON.
"""

import gzip
import hashlib
import json
import re
import sys
import time
from pathlib import Path
from urllib.parse import parse_qs, quote, unquote, urljoin, urlparse

import requests
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from logging_utils import get_logger

# ── Snakemake interface ───────────────────────────────────────
log = get_logger("download_tool_datasets", snakemake.log[0])
tool_sources_yaml = snakemake.input.tool_sources
out_manifest = snakemake.output.manifest
datasets_dir = Path(snakemake.params.datasets_dir)

# ── Constants ────────────────────────────────────────────────
REQUEST_TIMEOUT = 30  # seconds per HTTP request
RETRY_WAIT = 5  # seconds between retries
MAX_RETRIES = 3
GITHUB_API = "https://api.github.com"
OSF_API = "https://api.osf.io/v2"
SF_API = "https://sourceforge.net/projects/{project}/files/{path}/?format=json"

FASTA_EXTENSIONS = {".fa", ".fasta", ".faa", ".fa.gz", ".fasta.gz", ".faa.gz", ".zip"}


# ── HTTP helper ──────────────────────────────────────────────


def _get(url: str, **kwargs) -> requests.Response | None:
    """GET with retry and graceful failure. Returns None on any error."""
    headers = kwargs.pop("headers", {})
    headers.setdefault("User-Agent", "rnachallenge-pipeline/1.0")
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = requests.get(url, headers=headers, timeout=REQUEST_TIMEOUT, **kwargs)
            resp.raise_for_status()
            return resp
        except requests.RequestException as exc:
            log.warning(f"  Attempt {attempt}/{MAX_RETRIES} failed for {url}: {exc}")
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_WAIT)
    return None


def _save(content: bytes, dest: Path) -> bool:
    """Write bytes to dest; return True on success."""
    try:
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(content)
        return True
    except OSError as exc:
        log.warning(f"  Could not write {dest}: {exc}")
        return False


def _already_downloaded(dest: Path, url: str, tool: str) -> dict | None:
    """
    If dest exists and is non-empty, return a pre-built 'ok' result dict
    so the caller can skip the download.  Returns None if download is needed.
    """
    if dest.exists() and dest.stat().st_size > 0:
        log.info(f"  [{tool}] Already exists, skipping: {dest.name}")
        return {"tool": tool, "file": str(dest), "url": url, "status": "ok"}
    return None


def _safe_filename(url: str, fallback: str) -> str:
    """Derive a safe filename from a URL or use fallback."""
    parsed = urlparse(url)
    name = Path(parsed.path).name
    query_names = parse_qs(parsed.query).get("filename", [])
    if query_names:
        candidate = query_names[0].strip()
        if candidate:
            name = candidate
    # Fall back to query-string filenames when the path points to a generic
    # download script such as download.php or download_file.php.
    elif (
        not name
        or "." not in name
        or name.lower()
        in {
            "download",
            "download.php",
            "download_file",
            "download_file.php",
            "download.cgi",
        }
    ):
        qs = parsed.query
        m = re.search(r"filename=([^&]+)", qs)
        if m:
            name = unquote(m.group(1))
    return name if name and "." in name else fallback


def _dedupe_dest_name(dest: Path, url: str, used_names: set[str]) -> Path:
    """Ensure the destination filename is unique within the current download batch."""
    if dest.name not in used_names:
        used_names.add(dest.name)
        return dest

    suffix = hashlib.md5(url.encode()).hexdigest()[:8]
    candidate = dest.with_name(f"{dest.stem}_{suffix}{dest.suffix}")
    used_names.add(candidate.name)
    return candidate


def _is_fasta_like(name: str, allowed_extensions: list[str]) -> bool:
    """True if the filename looks like a sequence file."""
    name_lower = name.lower()
    for ext in allowed_extensions:
        if name_lower.endswith(ext.lower()):
            return True
    return False


def _extract_sourceforge_files_from_html(
    html_text: str,
    project: str,
    sf_path: str,
    allowed_extensions: list[str],
) -> list[str]:
    """
    Extract file paths from SourceForge HTML listings when the JSON endpoint
    serves an HTML page instead of API JSON.
    """
    results: list[str] = []
    seen: set[str] = set()
    norm_path = sf_path.strip("/")
    path_prefix = f"{norm_path}/" if norm_path else ""

    patterns = [
        r'href=["\']([^"\']+)["\']',
        r'data-url=["\']([^"\']+)["\']',
    ]

    for pattern in patterns:
        for raw_link in re.findall(pattern, html_text, flags=re.IGNORECASE):
            abs_link = urljoin("https://sourceforge.net", raw_link)
            parsed = urlparse(abs_link)
            path = parsed.path

            matched_root = None
            for root in (
                f"/project/{project}/files/",
                f"/projects/{project}/files/",
            ):
                if path.startswith(root):
                    matched_root = root
                    break

            if not matched_root:
                continue

            tail = path[len(matched_root) :].strip("/")
            if path_prefix and not tail.startswith(path_prefix):
                continue
            if path_prefix:
                tail = tail[len(path_prefix) :]

            if "/download" not in tail:
                continue

            rel_file = unquote(tail.split("/download", 1)[0].strip("/"))
            if not rel_file:
                continue
            if not _is_fasta_like(Path(rel_file).name, allowed_extensions):
                continue
            if rel_file in seen:
                continue
            seen.add(rel_file)
            results.append(rel_file)

    return results


def _extract_sourceforge_files_from_rss(
    rss_text: str,
    project: str,
    sf_path: str,
    allowed_extensions: list[str],
    include_all_files: bool = False,
) -> list[str]:
    """Extract SourceForge file paths from RSS item links."""
    results: list[str] = []
    seen: set[str] = set()
    norm_path = sf_path.strip("/")
    path_prefix = f"{norm_path}/" if norm_path else ""

    links = re.findall(r"<link>(.*?)</link>", rss_text, flags=re.IGNORECASE)
    for link in links:
        parsed = urlparse(link.strip())
        path = parsed.path

        matched_root = None
        for root in (
            f"/project/{project}/files/",
            f"/projects/{project}/files/",
        ):
            if path.startswith(root):
                matched_root = root
                break

        if not matched_root:
            continue

        tail = path[len(matched_root) :].strip("/")
        if path_prefix and not tail.startswith(path_prefix):
            continue
        if path_prefix:
            tail = tail[len(path_prefix) :]

        if "/download" not in tail:
            continue

        rel_file = unquote(tail.split("/download", 1)[0].strip("/"))
        if not rel_file:
            continue
        if not include_all_files and not _is_fasta_like(
            Path(rel_file).name, allowed_extensions
        ):
            continue
        if rel_file in seen:
            continue

        seen.add(rel_file)
        results.append(rel_file)

    return results


# ── Strategy implementations ─────────────────────────────────


def _download_github(tool: str, cfg: dict, out_dir: Path) -> list[dict]:
    """
    Enumerate files in a GitHub repo path via the Contents API and download them.
    Recursively handles subdirectories (up to depth 3).
    Supports either a single `path` or multiple `paths` entries.
    """
    owner = cfg["owner"]
    repo = cfg["repo"]
    branch = cfg.get("branch", "main")
    path = cfg.get("path", "")
    paths = cfg.get("paths", [])
    extensions = cfg.get("extensions", list(FASTA_EXTENSIONS))
    results = []

    # Backward compatible path handling:
    # - If `paths` is configured, download from each listed path.
    # - Otherwise use legacy single `path`.
    paths_to_process = [p for p in paths if p] if paths else [path]

    def _list_and_download(api_path: str, depth: int = 0):
        if depth > 3:
            return
        encoded_path = quote(api_path, safe="/")
        url = f"{GITHUB_API}/repos/{owner}/{repo}/contents/{encoded_path}?ref={branch}"
        resp = _get(url, headers={"Accept": "application/vnd.github+json"})
        if resp is None:
            log.warning(f"  [{tool}] GitHub API failed for path: {api_path}")
            return
        try:
            items = resp.json()
        except ValueError:
            log.warning(f"  [{tool}] GitHub API returned non-JSON for path: {api_path}")
            return
        if isinstance(items, dict) and "message" in items:
            log.warning(f"  [{tool}] GitHub API error: {items['message']}")
            return
        for item in items:
            if item.get("type") == "dir":
                _list_and_download(item["path"], depth + 1)
            elif item.get("type") == "file":
                name = item["name"]
                if not _is_fasta_like(name, extensions):
                    continue
                dl_url = item.get("download_url")
                if not dl_url:
                    continue
                dest = out_dir / name
                if cached := _already_downloaded(dest, dl_url, tool):
                    results.append(cached)
                    continue
                log.info(f"  [{tool}] Downloading {name} ...")
                resp2 = _get(dl_url)
                if resp2 and _save(resp2.content, dest):
                    results.append(
                        {"tool": tool, "file": str(dest), "url": dl_url, "status": "ok"}
                    )
                else:
                    results.append(
                        {
                            "tool": tool,
                            "file": str(dest),
                            "url": dl_url,
                            "status": "failed",
                        }
                    )

    for path_item in paths_to_process:
        _list_and_download(path_item)
    return results


def _download_direct_list(tool: str, cfg: dict, out_dir: Path) -> list[dict]:
    """Download each URL in the list directly."""
    results = []
    used_names: set[str] = set()
    for url in cfg.get("urls", []):
        name = _safe_filename(
            url, f"file_{hashlib.md5(url.encode()).hexdigest()[:8]}.fa"
        )
        dest = _dedupe_dest_name(out_dir / name, url, used_names)
        if cached := _already_downloaded(dest, url, tool):
            results.append(cached)
            continue
        log.info(f"  [{tool}] Downloading {name} from {url} ...")
        resp = _get(url)
        if resp and _save(resp.content, dest):
            results.append(
                {"tool": tool, "file": str(dest), "url": url, "status": "ok"}
            )
        else:
            results.append({"tool": tool, "file": "", "url": url, "status": "failed"})
    return results


def _download_sourceforge(tool: str, cfg: dict, out_dir: Path) -> list[dict]:
    """
    Use SourceForge JSON API to list files in one or more project directories,
    then download matching FASTA-like files.

    Config accepts either a single path string (sf_path) or a list of paths
    (sf_paths).  Both forms are normalised to a list internally.
    """
    project = cfg["project"]
    extensions = cfg.get("extensions", list(FASTA_EXTENSIONS))
    rss_include_all = bool(cfg.get("rss_include_all_files", False))
    results = []

    # Accept sf_path (string) or sf_paths (list); normalise to list
    raw = cfg.get("sf_paths") or cfg.get("sf_path", "")
    if isinstance(raw, str):
        paths = [raw.strip("/")]
    else:
        paths = [p.strip("/") for p in raw]

    for sf_path in paths:
        api_url = SF_API.format(project=project, path=sf_path)
        log.info(f"  [{tool}] Fetching SourceForge listing: {api_url}")
        resp = _get(api_url, headers={"Accept": "application/json"})
        if resp is None:
            log.warning(
                f"  [{tool}] SourceForge API unavailable for path '{sf_path}' — skip"
            )
            results.append(
                {
                    "tool": tool,
                    "file": "",
                    "url": api_url,
                    "status": "skipped",
                    "note": f"Manual download from https://sourceforge.net/projects/{project}/files/{sf_path}/",
                }
            )
            continue

        try:
            data = resp.json()
        except ValueError:
            fallback_source = ""
            preview = resp.text[:300].replace("\n", " ")
            log.warning(
                f"  [{tool}] SourceForge API returned non-JSON for path '{sf_path}' "
                f"(HTTP {resp.status_code}, Content-Type: {resp.headers.get('Content-Type', '?')}). "
                f"Response preview: {preview!r}"
            )
            fallback_files = _extract_sourceforge_files_from_html(
                resp.text,
                project,
                sf_path,
                extensions,
            )
            if fallback_files:
                fallback_source = "html"

            if not fallback_files:
                rss_url = (
                    f"https://sourceforge.net/projects/{project}/rss?path=/{sf_path}"
                )
                log.info(f"  [{tool}] Trying SourceForge RSS fallback: {rss_url}")
                rss_resp = _get(rss_url)
                if rss_resp is not None:
                    fallback_files = _extract_sourceforge_files_from_rss(
                        rss_resp.text,
                        project,
                        sf_path,
                        extensions,
                        include_all_files=rss_include_all,
                    )
                    if fallback_files:
                        fallback_source = "rss"

            if not fallback_files:
                results.append(
                    {
                        "tool": tool,
                        "file": "",
                        "url": api_url,
                        "status": "failed",
                        "note": (
                            f"SourceForge listing for '{sf_path}' is not JSON and no files "
                            "were discovered via HTML/RSS fallback"
                        ),
                    }
                )
                continue

            log.info(
                f"  [{tool}] SourceForge fallback discovered "
                f"{len(fallback_files)} candidate file(s)"
            )
            files = [
                {
                    "name": rel_file,
                    "from_sourceforge_fallback": True,
                    "sourceforge_fallback_source": fallback_source,
                }
                for rel_file in fallback_files
            ]
        else:
            files = data.get("files", [])

        if not files:
            log.warning(
                f"  [{tool}] SourceForge returned empty listing for path '{sf_path}'"
            )

        for item in files:
            name = item.get("name", "")
            from_fallback = bool(item.get("from_sourceforge_fallback", False))
            fallback_source = item.get("sourceforge_fallback_source", "")
            allow_non_fasta_from_rss = (
                from_fallback and fallback_source == "rss" and rss_include_all
            )

            if not allow_non_fasta_from_rss and not _is_fasta_like(name, extensions):
                continue
            encoded_name = quote(name, safe="/")
            dl_url = f"https://downloads.sourceforge.net/project/{project}/{sf_path}/{encoded_name}"
            dest = out_dir / Path(name).name
            if cached := _already_downloaded(dest, dl_url, tool):
                results.append(cached)
                continue
            log.info(f"  [{tool}] Downloading {name} ...")
            resp2 = _get(dl_url, allow_redirects=True)
            if resp2 and _save(resp2.content, dest):
                results.append(
                    {"tool": tool, "file": str(dest), "url": dl_url, "status": "ok"}
                )
            else:
                results.append(
                    {"tool": tool, "file": "", "url": dl_url, "status": "failed"}
                )

    return results


def _download_osf(tool: str, cfg: dict, out_dir: Path) -> list[dict]:
    """
    Use OSF API to enumerate files in a project's osfstorage and download
    any FASTA-like files.
    """
    node_id = cfg["node_id"]
    extensions = cfg.get("extensions", list(FASTA_EXTENSIONS))
    results = []
    next_url = f"{OSF_API}/nodes/{node_id}/files/osfstorage/"

    while next_url:
        resp = _get(next_url)
        if resp is None:
            log.warning(f"  [{tool}] OSF API unavailable")
            results.append(
                {
                    "tool": tool,
                    "file": "",
                    "url": next_url,
                    "status": "skipped",
                    "note": f"Manual download from https://osf.io/{node_id}/",
                }
            )
            break
        try:
            data = resp.json()
        except ValueError:
            log.warning(f"  [{tool}] OSF API returned non-JSON")
            break

        for item in data.get("data", []):
            attr = item.get("attributes", {})
            name = attr.get("name", "")
            kind = attr.get("kind", "")
            if kind == "folder":
                # Queue subfolder
                links = item.get("relationships", {}).get("files", {}).get("links", {})
                related = links.get("related", {}).get("href", "")
                if related:
                    # Recurse into subfolder (simple one-level expansion)
                    resp2 = _get(related)
                    if resp2:
                        try:
                            sub_data = resp2.json()
                            for sub_item in sub_data.get("data", []):
                                sub_attr = sub_item.get("attributes", {})
                                sub_name = sub_attr.get("name", "")
                                if not _is_fasta_like(sub_name, extensions):
                                    continue
                                sub_links = sub_item.get("links", {})
                                dl_url = sub_links.get("download", "")
                                if not dl_url:
                                    continue
                                dest = out_dir / sub_name
                                if cached := _already_downloaded(dest, dl_url, tool):
                                    results.append(cached)
                                    continue
                                log.info(f"  [{tool}] Downloading {sub_name} ...")
                                dl_resp = _get(dl_url)
                                if dl_resp and _save(dl_resp.content, dest):
                                    results.append(
                                        {
                                            "tool": tool,
                                            "file": str(dest),
                                            "url": dl_url,
                                            "status": "ok",
                                        }
                                    )
                                else:
                                    results.append(
                                        {
                                            "tool": tool,
                                            "file": "",
                                            "url": dl_url,
                                            "status": "failed",
                                        }
                                    )
                        except ValueError:
                            pass
                continue
            if not _is_fasta_like(name, extensions):
                continue
            links = item.get("links", {})
            dl_url = links.get("download", "")
            if not dl_url:
                continue
            dest = out_dir / name
            if cached := _already_downloaded(dest, dl_url, tool):
                results.append(cached)
                continue
            log.info(f"  [{tool}] Downloading {name} ...")
            dl_resp = _get(dl_url)
            if dl_resp and _save(dl_resp.content, dest):
                results.append(
                    {"tool": tool, "file": str(dest), "url": dl_url, "status": "ok"}
                )
            else:
                results.append(
                    {"tool": tool, "file": "", "url": dl_url, "status": "failed"}
                )

        # Follow pagination
        links = data.get("links", {})
        next_url = links.get("next")

    return results


def _download_readme_dir(tool: str, cfg: dict, out_dir: Path) -> list[dict]:
    """
    Download README files from a server directory, parse them for FASTA
    filenames, then download those.

    Fast-fail behaviour: if the first README probe returns a 404 (or any
    HTTP error), the entire server directory is considered gone and we
    immediately emit a manual_required result without probing any further
    URLs.  This avoids burning minutes retrying dozens of guessed filenames
    against a dead server.
    """
    base_url = cfg["base_url"].rstrip("/")
    subdirs = cfg.get("subdirs", [])
    readme_files = cfg.get("readme_files", [])
    extensions = cfg.get("extensions", list(FASTA_EXTENSIONS))
    results = []

    # ── Probe the first README to check whether the server is alive ──────
    if subdirs:
        probe_subdir = subdirs[0]
        probe_name = readme_files[0] if readme_files else f"README_{probe_subdir}"
        probe_url = f"{base_url}/{probe_subdir}/{probe_name}"
        log.info(f"  [{tool}] Probing server with: {probe_url}")
        probe_resp = _get(probe_url)
        if probe_resp is None:
            log.warning(
                f"  [{tool}] Server unreachable or returned an error — "
                "skipping all filename probing. MANUAL DOWNLOAD REQUIRED.\n"
                f"  Place .fa/.fasta files in {out_dir}"
            )
            return [
                {
                    "tool": tool,
                    "file": "",
                    "url": base_url,
                    "status": "manual_required",
                    "note": (
                        f"Server at {base_url} appears offline (404/connection error). "
                        f"Download manually and place .fa/.fasta files in {out_dir}"
                    ),
                }
            ]
    else:
        probe_resp = None

    # ── Server is alive — collect FASTA URLs from all READMEs ────────────
    all_fasta_urls = []

    for i, subdir in enumerate(subdirs):
        dir_url = f"{base_url}/{subdir}"
        readme_name = readme_files[i] if i < len(readme_files) else f"README_{subdir}"
        readme_url = f"{dir_url}/{readme_name}"

        # Reuse already-fetched response for the first subdir
        if i == 0 and probe_resp is not None:
            resp = probe_resp
        else:
            log.info(f"  [{tool}] Fetching README: {readme_url}")
            resp = _get(readme_url)

        if resp:
            for line in resp.text.splitlines():
                for ext in [".fa", ".fasta", ".gz", ".bed", ".gtf"]:
                    if ext in line.lower():
                        for tok in re.findall(
                            r"\S+" + re.escape(ext) + r"\S*", line, re.IGNORECASE
                        ):
                            tok = tok.strip(".,;:\"'()")
                            if "/" not in tok:
                                all_fasta_urls.append(f"{dir_url}/{tok}")
        else:
            log.warning(
                f"  [{tool}] README not accessible at {readme_url} — skipping subdir"
            )

    # ── Download discovered files ─────────────────────────────────────────
    all_fasta_urls = list(dict.fromkeys(all_fasta_urls))  # deduplicate, preserve order

    for url in all_fasta_urls:
        name = Path(urlparse(url).path).name
        if not _is_fasta_like(name, extensions):
            continue
        dest = out_dir / name
        if cached := _already_downloaded(dest, url, tool):
            results.append(cached)
            continue
        log.info(f"  [{tool}] Downloading {name} ...")
        resp = _get(url)
        if resp and _save(resp.content, dest):
            results.append(
                {"tool": tool, "file": str(dest), "url": url, "status": "ok"}
            )
        else:
            results.append({"tool": tool, "file": "", "url": url, "status": "failed"})

    if not any(r["status"] == "ok" for r in results):
        log.warning(
            f"  [{tool}] No files downloaded. Place .fa/.fasta files manually in {out_dir}"
        )
        results.append(
            {
                "tool": tool,
                "file": "",
                "url": base_url,
                "status": "skipped",
                "note": f"Manual download needed — place .fa/.fasta files in {out_dir}",
            }
        )

    return results


def _download_web(tool: str, cfg: dict, out_dir: Path) -> list[dict]:
    """
    Attempt a single URL. Most 'web' strategy tools are portal pages that
    require human navigation — log the manual download instruction.
    """
    url = cfg["url"]
    note = cfg.get("note", "")
    log.warning(
        f"  [{tool}] Strategy 'web': attempting {url}\n" f"  If this fails, {note}"
    )
    resp = _get(url)
    if resp is None:
        log.warning(f"  [{tool}] Could not access {url}. MANUAL DOWNLOAD REQUIRED.")
        manual_note = (
            note
            if note
            else (f"Please download manually and place .fa/.fasta files in {out_dir}")
        )
        return [
            {
                "tool": tool,
                "file": "",
                "url": url,
                "status": "manual_required",
                "note": manual_note,
            }
        ]

    # If we got a response, check if it's a FASTA (Content-Type or content check)
    content_type = resp.headers.get("Content-Type", "")
    if "text/html" in content_type:
        # It's a webpage — not directly downloadable
        log.warning(
            f"  [{tool}] {url} returned HTML (portal page). MANUAL DOWNLOAD REQUIRED.\n"
            f"  {note}"
        )
        return [
            {
                "tool": tool,
                "file": "",
                "url": url,
                "status": "manual_required",
                "note": note,
            }
        ]

    # Looks like a direct file
    name = _safe_filename(url, f"{tool}_download.fa")
    dest = out_dir / name
    if _save(resp.content, dest):
        return [{"tool": tool, "file": str(dest), "url": url, "status": "ok"}]
    return [{"tool": tool, "file": "", "url": url, "status": "failed"}]


# ── Dispatch ────────────────────────────────────────────────

STRATEGY_FUNCS = {
    "github": _download_github,
    "direct_list": _download_direct_list,
    "sourceforge": _download_sourceforge,
    "osf": _download_osf,
    "readme_dir": _download_readme_dir,
    "web": _download_web,
}


# ── Main ────────────────────────────────────────────────────

log.info("Stage 0: Downloading tool dataset FASTAs")

with open(tool_sources_yaml) as fh:
    sources = yaml.safe_load(fh)

tools_cfg = sources.get("tools", {})
manifest: list[dict] = []
manual_required: list[str] = []

for tool_name, cfg in tools_cfg.items():
    strategy = cfg.get("strategy", "web")
    out_dir = datasets_dir / tool_name
    out_dir.mkdir(parents=True, exist_ok=True)

    log.info(f"[{tool_name}] strategy={strategy}")

    func = STRATEGY_FUNCS.get(strategy)
    if func is None:
        log.warning(f"[{tool_name}] Unknown strategy '{strategy}' — skipping")
        manifest.append(
            {"tool": tool_name, "file": "", "url": "", "status": "unknown_strategy"}
        )
        continue

    try:
        results = func(tool_name, cfg, out_dir)
    except Exception as exc:
        log.error(f"[{tool_name}] Unexpected error: {exc}", exc_info=True)
        results = [
            {
                "tool": tool_name,
                "file": "",
                "url": "",
                "status": "error",
                "note": str(exc),
            }
        ]

    manifest.extend(results)

    ok_count = sum(1 for r in results if r["status"] == "ok")

    # Reconcile with local files: any non-empty tool folder is considered
    # usable for downstream stages, regardless of remote download status.
    local_files = sorted([p for p in out_dir.rglob("*") if p.is_file()])
    existing_ok_files = {
        Path(r["file"]).resolve()
        for r in results
        if r.get("status") == "ok" and r.get("file")
    }
    added_local_ok = 0
    if local_files:
        for fp in local_files:
            resolved = fp.resolve()
            if resolved in existing_ok_files:
                continue
            manifest.append(
                {
                    "tool": tool_name,
                    "file": str(fp),
                    "url": "",
                    "status": "ok",
                    "note": "local_file_detected",
                }
            )
            added_local_ok += 1

        if added_local_ok:
            log.info(
                f"[{tool_name}] Added {added_local_ok} local file(s) as status=ok "
                "for downstream processing"
            )

        # Local files present means manual intervention has already been done.
        if tool_name in manual_required:
            manual_required = [t for t in manual_required if t != tool_name]

    ok_count += added_local_ok
    log.info(f"[{tool_name}] {ok_count}/{len(results)} files downloaded OK")

    if len(results) == 0:
        log.warning(
            f"[{tool_name}] No files were found at the configured location "
            f"(strategy='{strategy}'). This is likely a configuration problem — "
            f"check the path, extensions, or repository name in config/tool_sources.yaml "
            f"for '{tool_name}'."
        )
        manifest.append(
            {
                "tool": tool_name,
                "file": "",
                "url": "",
                "status": "config_error",
                "note": (
                    f"No files matched for strategy='{strategy}'. "
                    "Review path/extensions in config/tool_sources.yaml."
                ),
            }
        )
        manual_required.append(tool_name)

    if any(r["status"] == "manual_required" for r in results):
        if not local_files:
            manual_required.append(tool_name)

# ── Write manifest ───────────────────────────────────────────
Path(out_manifest).parent.mkdir(parents=True, exist_ok=True)
with open(out_manifest, "w") as fh:
    json.dump(manifest, fh, indent=2)

log.info("=" * 60)
log.info(f"Download manifest written to: {out_manifest}")

ok_total = sum(1 for r in manifest if r["status"] == "ok")
fail_total = sum(1 for r in manifest if r["status"] == "failed")
skip_total = sum(1 for r in manifest if r["status"] in ("skipped", "manual_required"))
config_total = sum(1 for r in manifest if r["status"] == "config_error")

log.info(f"  Downloaded OK      : {ok_total}")
log.info(f"  Failed             : {fail_total}")
if config_total:
    log.warning(f"  Config errors      : {config_total}  ← check tool_sources.yaml")
log.info(f"  Manual required    : {skip_total}")

if manual_required:
    log.warning("=" * 60)
    log.warning("MANUAL DOWNLOAD REQUIRED for the following tools:")
    for t in manual_required:
        note = next(
            (r.get("note", "") for r in manifest if r["tool"] == t and r.get("note")),
            "",
        )
        log.warning(f"  {t}: {note}")
    log.warning(
        "Place the FASTA files in resources/tool_datasets/<ToolName>/ "
        "and rerun the pipeline."
    )

log.info("Stage 0 complete.")
