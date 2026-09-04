#!/bin/bash

# ==========================================================
# Module 4B
# BioSample -> SRA RunInfo retrieval
# ==========================================================

set -u

INPUT="outputs/04_ast_biosamples.txt"
RAW_OUTPUT="temp/04B_runinfo_raw.csv"
LOG="logs/04B_retrieval.log"
FAILED="logs/04B_failed_biosamples.txt"

mkdir -p outputs logs temp

if [ ! -f "$INPUT" ]; then
    echo "ERROR: Missing input file: $INPUT"
    exit 1
fi

rm -f "$RAW_OUTPUT" "$LOG" "$FAILED"

TOTAL=$(grep -cve '^[[:space:]]*$' "$INPUT")
COUNT=0
SUCCESS=0
NO_RUNS=0

echo "Starting Module 4B..." | tee "$LOG"
echo "BioSamples to process: $TOTAL" | tee -a "$LOG"

# Read from file descriptor 3.
# This prevents EDirect from consuming the BioSample input file.
while IFS= read -r BIOSAMPLE <&3 || [ -n "$BIOSAMPLE" ]
do
    BIOSAMPLE=$(echo "$BIOSAMPLE" | tr -d '\r' | xargs)

    if [ -z "$BIOSAMPLE" ]; then
        continue
    fi

    COUNT=$((COUNT + 1))

    echo "[$COUNT/$TOTAL] $BIOSAMPLE" | tee -a "$LOG"

    TEMP_FILE="temp/runinfo_${BIOSAMPLE}.csv"

    rm -f "$TEMP_FILE"

    esearch \
        -db biosample \
        -query "${BIOSAMPLE}[Accession]" \
    | elink \
        -target sra \
    | efetch \
        -format runinfo \
        > "$TEMP_FILE" 2>> "$LOG"

    if [ -s "$TEMP_FILE" ] && [ "$(wc -l < "$TEMP_FILE")" -gt 1 ]; then

        if [ ! -s "$RAW_OUTPUT" ]; then
            cat "$TEMP_FILE" > "$RAW_OUTPUT"
        else
            tail -n +2 "$TEMP_FILE" >> "$RAW_OUTPUT"
        fi

        SUCCESS=$((SUCCESS + 1))

    else
        echo "$BIOSAMPLE" >> "$FAILED"
        NO_RUNS=$((NO_RUNS + 1))
        echo "  No SRA RunInfo found" | tee -a "$LOG"
    fi

    rm -f "$TEMP_FILE"

    sleep 0.34

done 3< "$INPUT"

echo "" | tee -a "$LOG"
echo "==============================================" | tee -a "$LOG"
echo "Module 4B retrieval completed" | tee -a "$LOG"
echo "BioSamples processed: $COUNT" | tee -a "$LOG"
echo "BioSamples with SRA RunInfo: $SUCCESS" | tee -a "$LOG"
echo "BioSamples without SRA RunInfo: $NO_RUNS" | tee -a "$LOG"
echo "Raw RunInfo file: $RAW_OUTPUT" | tee -a "$LOG"
echo "Failed BioSamples: $FAILED" | tee -a "$LOG"
echo "==============================================" | tee -a "$LOG"


