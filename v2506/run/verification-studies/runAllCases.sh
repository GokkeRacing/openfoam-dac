#!/bin/bash
# Script to run all OpenFOAM verification cases in parallel

# Define the range of cases
START_CASE=0
END_CASE=24

echo "Starting parallel simulation of all cases (case_0000 to case_0028)..."
echo "================================================"

# Function to run a single case
run_case() {
    local CASE_DIR=$1
    local CASE_NUM=$2
    
    if [ ! -d "$CASE_DIR" ]; then
        echo "Warning: Directory $CASE_DIR does not exist. Skipping..."
        return 1
    fi
    
    echo "[$CASE_NUM] Starting $CASE_DIR..."
    
    cd "$CASE_DIR"
    
    # Run Allpre if it exists
    if [ -f "Allpre" ]; then
        ./Allpre > Allpre.log 2>&1
    else
        echo "[$CASE_NUM] Warning: Allpre not found in $CASE_DIR"
    fi
    
    # Run Allrun-serial if it exists
    if [ -f "Allrun-serial" ]; then
        ./Allrun-serial > Allrun.log 2>&1
    else
        echo "[$CASE_NUM] Warning: Allrun-serial not found in $CASE_DIR"
    fi
    
    cd ..
    
    echo "[$CASE_NUM] Completed $CASE_DIR"
}

# Export the function so it's available to subshells
export -f run_case

# Loop through all cases and launch them in parallel
for i in $(seq -f "%04g" $START_CASE $END_CASE); do
    CASE_DIR="case_$i"
    run_case "$CASE_DIR" "$i" &
done

# Wait for all background processes to complete
echo ""
echo "Waiting for all cases to complete..."
wait

echo ""
echo "================================================"
echo "All cases completed!"