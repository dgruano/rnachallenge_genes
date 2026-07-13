"""Branch tests for download_phytozome_gtf: a stale manifest ``status`` must
not gate the live JGI resolution (else PURGED entries can never self-heal on
rerun). No network — the JGI lookup functions are monkeypatched.
"""

import json
import sys
import types
from pathlib import Path

import pytest

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import download_phytozome_gtf as dl
from snakemake.exceptions import WorkflowError


def _snakemake(tmp_path, entry, config=None):
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"vitis_vinifera": entry}))
    return types.SimpleNamespace(
        wildcards=types.SimpleNamespace(species="vitis_vinifera"),
        input=types.SimpleNamespace(manifest=str(manifest)),
        output=[str(tmp_path / "out.gff3.gz")],
        log=[str(tmp_path / "out.log")],
        config=config or {},
    )


def test_stale_purged_with_genome_id_fires_live_restore(tmp_path, monkeypatch):
    """Stale PURGED + a genome_id → reach the live path and fire a restore."""
    calls = {}
    monkeypatch.setattr(
        dl,
        "resolve_annotation",
        lambda gid, tok, prefer_name=None: {"_id": "abc", "file_status": "PURGED"},
    )
    monkeypatch.setattr(
        dl,
        "request_restore",
        lambda ids, tok: calls.setdefault("restore", ids) or "queued",
    )
    monkeypatch.setattr(dl, "load_token", lambda: "tok")

    sm = _snakemake(
        tmp_path,
        {"status": "PURGED"},
        config={"phytozome_gtf_sources": {"vitis_vinifera": {"genome_id": 145}}},
    )
    with pytest.raises(WorkflowError, match="cold storage"):
        dl.main(sm)
    assert calls["restore"] == "abc"  # live restore actually fired


def test_missing_manifest_entry_falls_through_to_config_genome_id(
    tmp_path, monkeypatch
):
    """Species absent from manifest but present in config → resolve via config
    genome_id instead of hard-failing (the manifest only pins portal_file_name)."""
    seen = {}
    monkeypatch.setattr(dl, "load_token", lambda: "tok")

    def fake_resolve(gid, tok, prefer_name=None):
        seen["gid"] = gid
        return {"_id": "x", "file_status": "RESTORED", "download_url": "http://x"}

    monkeypatch.setattr(dl, "resolve_annotation", fake_resolve)

    def no_network(*a, **k):
        raise RuntimeError("stub: reached download")

    monkeypatch.setattr(dl.urllib.request, "urlopen", no_network)

    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"other_species": {"genome_id": 1}}))
    sm = types.SimpleNamespace(
        wildcards=types.SimpleNamespace(species="solanum_tuberosum"),
        input=types.SimpleNamespace(manifest=str(manifest)),
        output=[str(tmp_path / "out.gff3.gz")],
        log=[str(tmp_path / "out.log")],
        config={"phytozome_gtf_sources": {"solanum_tuberosum": {"genome_id": 206}}},
    )
    # urlopen is stubbed to fail; we only assert the config genome_id was reached.
    with pytest.raises(Exception):
        dl.main(sm)
    assert seen["gid"] == 206


def test_stale_purged_without_genome_id_reports_clearly(tmp_path, monkeypatch):
    """Stale PURGED + no genome_id → clean message, no restore attempted."""
    monkeypatch.setattr(dl, "load_token", lambda: "tok")
    monkeypatch.setattr(
        dl,
        "request_restore",
        lambda *a, **k: pytest.fail("should not fire restore without genome_id"),
    )
    sm = _snakemake(tmp_path, {"status": "PURGED"}, config={})
    with pytest.raises(WorkflowError, match="no genome_id"):
        dl.main(sm)
