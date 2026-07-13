#!/bin/bash
while true; do
  clear
  echo "=== NCBI Resolver Progress Monitor ==="
  echo "Time: $(date '+%H:%M:%S')"
  echo ""

  # Check if job is running
  if squeue -u $USER | grep -q 81920; then
    echo "Status: RUNNING (SLURM job 81920)"
  else
    echo "Status: COMPLETED or FAILED"
  fi
  echo ""

  # Count batches
  BATCH=$(grep "NCBI batch" logs/resolve_ids.log 2>/dev/null | tail -1 | grep -oP 'batch \K\d+')
  echo "Current batch: $BATCH / ~185"
  echo ""

  # Count resolution methods
  echo "Resolution Statistics:"
  EFETCH=$(grep -c "resolved via efetch fallback" logs/resolve_ids.log 2>/dev/null || echo 0)
  BASE=$(grep -c "resolved via base ID" logs/resolve_ids.log 2>/dev/null || echo 0)
  FAILED=$(grep -c "could not resolve.*even via efetch" logs/resolve_ids.log 2>/dev/null || echo 0)

  echo "  efetch fallback: $EFETCH"
  echo "  base ID search:  $BASE"
  echo "  failed:          $FAILED"
  echo ""

  # Recent log entries
  echo "Recent activity:"
  tail -5 logs/resolve_ids.log 2>/dev/null | sed 's/^/  /'

  sleep 15
done
