#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --time=4-00:00:00
#SBATCH --mem=64000mb
#SBATCH --job-name=validation5.1
#SBATCH --error=job.%A.err
#SBATCH --output=job.%A.out

# Perform the operation
while IFS= read -r line
do
  if [ "$line" != "" ]; then
    python3 pp2_predictor.py "$line"
  fi
done < <(tail -n +2 input_database.txt)
