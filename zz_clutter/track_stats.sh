#!/bin/bash
# Track resolution statistics during run

LOG_FILE="/mnt/cbib/LNClassifier/paper/rnachallenge_genes/logs/resolve_external_ids.log"
STATS_FILE="/mnt/cbib/LNClassifier/paper/rnachallenge_genes/logs/resolution_stats.txt"

echo "Timestamp,IDs_Processed,Resolved,Ambiguous,Unresolved,Resolution_Rate,Ambiguity_Rate" > "$STATS_FILE"

while true; do
    # Get latest progress line
    progress=$(tail -100 "$LOG_FILE" | grep "Progress:" | tail -1)

    if [ -n "$progress" ]; then
        # Extract numbers using regex
        ids_processed=$(echo "$progress" | grep -oP '\d+/\d+' | cut -d'/' -f1)
        resolved=$(echo "$progress" | grep -oP 'Resolved: \K\d+')
        ambiguous=$(echo "$progress" | grep -oP 'Ambiguous: \K\d+')
        unresolved=$(echo "$progress" | grep -oP 'Unresolved: \K\d+')

        if [ -n "$ids_processed" ] && [ -n "$resolved" ]; then
            res_rate=$(awk "BEGIN {printf \"%.1f\", 100*$resolved/$ids_processed}")
            amb_rate=$(awk "BEGIN {printf \"%.1f\", 100*$ambiguous/$ids_processed}")

            timestamp=$(date '+%Y-%m-%d %H:%M:%S')
            echo "$timestamp,$ids_processed,$resolved,$ambiguous,$unresolved,$res_rate,$amb_rate" >> "$STATS_FILE"
        fi
    fi

    # Check if complete
    if grep -q "Resolution complete" "$LOG_FILE" 2>/dev/null; then
        echo "Resolution completed at $(date)" >> "$STATS_FILE"
        break
    fi

    # Check if process still running
    if ! pgrep -f "resolve_external_ids.py" > /dev/null; then
        echo "Process not found at $(date)" >> "$STATS_FILE"
        break
    fi

    sleep 120  # Check every 2 minutes
done
