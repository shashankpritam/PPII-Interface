#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --time=4-00:00:00
#SBATCH --mem=64000mb
#SBATCH --job-name=validation5.1
#SBATCH --error=job.%A.err
#SBATCH --output=job.%A.out

# Create a CSV file and write the headers
echo "PDB_ID,Model,Chain,TRP,NBR,Template_PPII,Score" > output.csv
echo "PDB_ID,Model,Chain,TRP,NBR,Template_PPII,Score" > binding_site_output.csv

# Perform the operation
while IFS= read -r line
do
  if [ "$line" != "" ]; then
    # Run the Python script and save the output in a variable
    output=$(python3 pp2_predictor.py "$line")
    
    # Extract needed information using regex and echo to the csv file
    if [[ "$output" =~ (For query structure )([^,]+)(, predicted binding site details are - Model = )([^,]+)(, Chain = )([^,]+)(, TRP = )([^,]+)(, NBR = )([^,]+)(Template PPII is \[')([^']+)('\], with a Score of )([^$]+) ]]; then
      PDB_ID=${BASH_REMATCH[2]}
      Model=${BASH_REMATCH[4]}
      Chain=${BASH_REMATCH[6]}
      TRP=${BASH_REMATCH[8]}
      NBR=${BASH_REMATCH[10]}
      Template_PPII=${BASH_REMATCH[12]}
      Score=${BASH_REMATCH[14]}
      echo "$PDB_ID,$Model,$Chain,$TRP,$NBR,$Template_PPII,$Score" >> output.csv
      echo "$PDB_ID,$Model,$Chain,$TRP,$NBR,$Template_PPII,$Score" >> binding_site_output.csv
    fi
  fi
done < <(tail -n +2 input_database.txt)
