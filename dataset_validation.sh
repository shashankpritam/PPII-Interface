#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --time=8-00:00:00
#SBATCH --mem=64000mb
#SBATCH --job-name=pp2_validation
#SBATCH --error=job.%A.err
#SBATCH --output=job.%A.out

# Defining the CSV headers
CSV_HEADERS="PDB_ID,Model,Chain,TRP,NBR,Template_PPII,Score"

# Output files
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
  # Skip empty lines
  if [ "$pdb_id" != "" ]; then
    # Run the Python script and save the output in a variable
    output=$(python3 pp2_predictor.py "$pdb_id")
    
    # If the output matches the expected pattern, extract information and write to the CSV files
    if [[ "$output" =~ (For query structure )([^,]+)(, predicted binding site details are - Model = )([^,]+)(, Chain = )([^,]+)(, TRP = )([^,]+)(, NBR = )([^,]+)(Template PPII is \\\[\\\')([^\'\\]]+)(\\\', with a Score of )([^$]+) ]]; then
      PDB_ID=${BASH_REMATCH[2]}
      Model=${BASH_REMATCH[4]}
      Chain=${BASH_REMATCH[6]}
      TRP=${BASH_REMATCH[8]}
      NBR=${BASH_REMATCH[10]}
      Template_PPII=${BASH_REMATCH[12]}
      Score=${BASH_REMATCH[14]}
      
      # Construct the CSV line
      CSV_LINE="$PDB_ID,$Model,$Chain,$TRP,$NBR,$Template_PPII,$Score"
      
      # Write the CSV line to the output files
      echo $CSV_LINE >> $OUTPUT_FILE
      echo $CSV_LINE >> $BINDING_SITE_FILE
    fi
  fi
done <<< "$TAIL_DATA_FILE"
