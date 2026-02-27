#!/bin/bash
# Autonomous monitoring script for external ID resolver
# Checks progress every 5 minutes and stops if errors detected

LOG_FILE="/mnt/cbib/LNClassifier/paper/rnachallenge_genes/logs/resolve_external_ids.log"
MONITOR_LOG="/mnt/cbib/LNClassifier/paper/rnachallenge_genes/logs/monitor.log"
SNAKEMAKE_LOG="/mnt/cbib/LNClassifier/paper/rnachallenge_genes/logs/manual_run.log"
CHECK_INTERVAL=300  # 5 minutes

echo "=== Monitoring started at $(date) ===" >> "$MONITOR_LOG"

last_progress=""
stall_count=0
MAX_STALLS=3  # Stop if no progress for 15 minutes (3 × 5min checks)

while true; do
    sleep $CHECK_INTERVAL
    
    echo "--- Check at $(date) ---" >> "$MONITOR_LOG"
    
    # Check if Snakemake process is still running
    if ! pgrep -f "snakemake.*resolve_external_ids" > /dev/null; then
        echo "Snakemake process not found - job may have completed or crashed" >> "$MONITOR_LOG"
        
        # Check if output files exist
        if [ -f "results/external_resolved.tsv" ] && \
           [ -f "results/external_ambiguous.tsv" ] && \
           [ -f "results/external_unresolved.tsv" ]; then
            echo "SUCCESS: All output files created!" >> "$MONITOR_LOG"
            tail -10 "$LOG_FILE" >> "$MONITOR_LOG"
            echo "=== Monitoring completed successfully at $(date) ===" >> "$MONITOR_LOG"
            exit 0
        else
            echo "ERROR: Output files missing - job may have failed" >> "$MONITOR_LOG"
            echo "Last 20 lines of resolver log:" >> "$MONITOR_LOG"
            tail -20 "$LOG_FILE" >> "$MONITOR_LOG"
            echo "Last 30 lines of Snakemake log:" >> "$MONITOR_LOG"
            tail -30 "$SNAKEMAKE_LOG" >> "$MONITOR_LOG"
            exit 1
        fi
    fi
    
    # Get latest progress
    current_progress=$(tail -5 "$LOG_FILE" | grep "Progress:" | tail -1)
    
    if [ -z "$current_progress" ]; then
        echo "WARNING: No progress line found in log" >> "$MONITOR_LOG"
    else
        echo "Current: $current_progress" >> "$MONITOR_LOG"
        
        # Check for stall (same progress as last check)
        if [ "$current_progress" == "$last_progress" ]; then
            ((stall_count++))
            echo "WARNING: No progress change detected (stall count: $stall_count/$MAX_STALLS)" >> "$MONITOR_LOG"
            
            if [ $stall_count -ge $MAX_STALLS ]; then
                echo "ERROR: Process appears stalled - no progress for $((MAX_STALLS * CHECK_INTERVAL / 60)) minutes" >> "$MONITOR_LOG"
                echo "Killing stalled processes..." >> "$MONITOR_LOG"
                pkill -f "snakemake.*resolve_external_ids"
                pkill -f "resolve_external_ids.py"
                echo "Last 30 lines of log:" >> "$MONITOR_LOG"
                tail -30 "$LOG_FILE" >> "$MONITOR_LOG"
                exit 1
            fi
        else
            stall_count=0
        fi
        
        last_progress="$current_progress"
    fi
    
    # Check for ERROR or WARNING in recent log lines
    if tail -20 "$LOG_FILE" | grep -iE "ERROR|Exception|Traceback" > /dev/null; then
        echo "ERROR detected in log file:" >> "$MONITOR_LOG"
        tail -30 "$LOG_FILE" >> "$MONITOR_LOG"
        echo "Stopping processes..." >> "$MONITOR_LOG"
        pkill -f "snakemake.*resolve_external_ids"
        pkill -f "resolve_external_ids.py"
        exit 1
    fi
    
    # Log warnings but don't stop (REST API warnings are expected)
    warning_count=$(tail -100 "$LOG_FILE" | grep -c "WARNING")
    if [ "$warning_count" -gt 0 ]; then
        echo "Found $warning_count warnings in last 100 log lines (expected for failed REST lookups)" >> "$MONITOR_LOG"
    fi
done
