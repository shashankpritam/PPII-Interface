#!/bin/bash

# SLURM job settings
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --time=8-00:00:00
#SBATCH --mem=64000mb
#SBATCH --job-name=pp2_validation
#SBATCH --error=job.%A.err
#SBATCH --output=job.%A.out

# Define the CSV headers and output files
CSV_HEADERS="PDB_ID,Model,Chain,TRP,NBR,Template_PPII,Score"
OUTPUT_FILE="output.csv"
BINDING_SITE_FILE="binding_site_output.csv"

# Initialize the CSV files with headers
echo $CSV_HEADERS > $OUTPUT_FILE
echo $CSV_HEADERS > $BINDING_SITE_FILE

# Exclude the first line of the input file (header line)
DATA_FILE="input_database.txt"
TAIL_DATA_FILE=$(tail -n +2 $DATA_FILE)

# Loop through each line of the data file
while IFS= read -r pdb_id
do
  # Only proceed for non-empty lines
  if [ "$pdb_id" != "" ]; then
    # Log the start time
    START_TIME=$(date +%s)

    # Run the Python script and capture its output
    output=$(python3 pp2_predictor.py "$pdb_id")

    # If the output matches the expected format, extract the details
    if [[ "$output" =~ (For query structure )([^,]+)(, predicted binding site details are - Model = )([^,]+)(, Chain = )([^,]+)(, TRP = )([^,]+)(, NBR = )([^,]+)(Template PPII is \\\[\\\')([^\'\\]]+)(\\\', with a Score of )([^$]+) ]]; then
      # Extract data using regex group matches
      PDB_ID=${BASH_REMATCH[2]}
      Model=${BASH_REMATCH[4]}
      Chain=${BASH_REMATCH[6]}
      TRP=${BASH_REMATCH[8]}
      NBR=${BASH_REMATCH[10]}
      Template_PPII=${BASH_REMATCH[12]}
      Score=${BASH_REMATCH[14]}
      
      # Construct the CSV line
      CSV_LINE="$PDB_ID,$Model,$Chain,$TRP,$NBR,$Template_PPII,$Score"
      
      # Append the line to both output files
      echo $CSV_LINE | tee -a $OUTPUT_FILE $BINDING_SITE_FILE
    fi

    # Calculate and log the elapsed time
    END_TIME=$(date +%s)
    ELAPSED_TIME=$((END_TIME - START_TIME))
    echo "Processing PDB ID $pdb_id took $ELAPSED_TIME seconds."
  fi
done <<< "$TAIL_DATA_FILE"
