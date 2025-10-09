#!/bin/bash
# -------------------------------
# run script

# just ./run.sh to run

# if shows error permission denied then first run
# chmod +x run.sh

# -------------------------------

set -e

# Ask if user wants a log file
while true; do
    read -p "Do you want to save a log file? (Y/N) " SAVE_LOG
    SAVE_LOG=$(echo "$SAVE_LOG" | tr '[:upper:]' '[:lower:]')
    if [[ "$SAVE_LOG" == "y" || "$SAVE_LOG" == "n" ]]; then
        break
    else
        echo "Invalid input. Enter Y or N."
    fi
done

if [[ "$SAVE_LOG" == "y" ]]; then
    LOG_FILE="log.run"
    LOG_CMD="| tee -a $LOG_FILE"
else
    LOG_CMD=""
fi

# Clean the case first
echo "Cleaning previous simulation data..."
foamCleanTutorials | tee -a $LOG_FILE  # cleans time directories, processor directories, logs, etc.

# Backup initial conditions
if [ -d "0.orig" ]; then
    echo "Copying initial conditions from 0.orig..."
    cp -r 0.orig 0
else
    echo "Error: 0.orig directory not found!"
    exit 1
fi

# Ask user for run mode
while true; do
    read -p "Run mode (serial/parallel)? " RUN_MODE
    RUN_MODE=$(echo "$RUN_MODE" | tr '[:upper:]' '[:lower:]')  # convert to lowercase
    if [[ "$RUN_MODE" == "serial" || "$RUN_MODE" == "parallel" ]]; then
        break
    else
        echo "Invalid input. Please enter 'serial' or 'parallel'."
    fi
done

# Run blockMesh
echo "Running blockMesh..."
eval "blockMesh $LOG_CMD"

# Run setFields
echo "Running setFields..."
eval "setFields $LOG_CMD"

# Run according to mode
if [ "$RUN_MODE" == "serial" ]; then
    echo "Running porousCapillaryFoam in SERIAL..."
    eval "porousCapillaryFoam $LOG_CMD"

elif [ "$RUN_MODE" == "parallel" ]; then
    # Ask for number of cores
    while true; do
        read -p "Enter number of cores (e.g., 2,3,4...): " NUM_PROCS
        if [[ "$NUM_PROCS" =~ ^[0-9]+$ && "$NUM_PROCS" -ge 2 ]]; then
            break
        else
            echo "Invalid input. Please enter an integer >= 2."
        fi
    done

    echo "Preparing decomposeParDict for $NUM_PROCS processors..."

    # Create hierarchical decomposeParDict
    cat > system/decomposeParDict <<EOL
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2306                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}

numberOfSubdomains $NUM_PROCS;

method  hierarchical;

coeffs
{
    n               (1 $NUM_PROCS 1);
}

EOL

    # Decompose the domain
    echo "Decomposing domain..."
    eval "decomposePar -force $LOG_CMD"

    # Run solver in parallel
    echo "Running porousCapillaryFoam in PARALLEL..."
    eval "mpirun -np $NUM_PROCS porousCapillaryFoam -parallel $LOG_CMD"

    # Reconstruct domain
    echo "Reconstructing domain..."
    eval "reconstructPar -latestTime $LOG_CMD"
fi

echo "Simulation completed successfully!"

