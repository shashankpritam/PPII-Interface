import os
import sys
import subprocess
import re
import csv

# Regular expression pattern for extracting the data
pattern = r"For query structure (.*?), predicted binding site details are - Model = (.*?), Chain = (.*?), TRP = (.*?), NBR = (.*?). Template PPII is \['(.*?)'\], with a Score of (.*?)$"

def run_pp2_predictor(pdb):
    # Construct the command
    command = ["python3", "pp2_predictor.py", pdb.strip()]

    try:
        # Execute the command and capture the output
        process = subprocess.run(command, check=True, text=True, capture_output=True)
        return process.stdout
    except subprocess.CalledProcessError as e:
        print(f"pp2_predictor.py failed for PDB ID {pdb.strip()}. Error: {e}", file=sys.stderr)
        return None

def main():
    # Check if pdb ids file is provided
    if len(sys.argv) < 2:
        print("Usage: python3 script.py input_database.txt", file=sys.stderr)
        sys.exit(1)

    pdb_ids_file = sys.argv[1]
    
    if not os.path.isfile(pdb_ids_file):
        print(f"File {pdb_ids_file} does not exist.", file=sys.stderr)
        sys.exit(1)

    # Open the CSV files for writing
    with open(pdb_ids_file, 'r') as file, open('output.csv', 'w', newline='') as output_csvfile, open('binding_site_output.csv', 'w', newline='') as binding_site_csvfile:
        output_writer = csv.writer(output_csvfile)
        binding_site_writer = csv.writer(binding_site_csvfile)

        # Write the headers to the CSV files
        headers = ["PDB_ID", "Model", "Chain", "TRP", "NBR", "Template_PPII", "Score"]
        output_writer.writerow(headers)
        binding_site_writer.writerow(headers)

        lines = file.readlines()[1:]  # Skip header line
        for pdb in lines:
            if pdb.strip():  # Ignore empty lines
                output = run_pp2_predictor(pdb)
                if output is not None:
                    # Extract the needed information using regex
                    match = re.search(pattern, output)
                    if match:
                        data = list(match.groups())
                        # Write the data to the CSV files
                        output_writer.writerow(data)
                        binding_site_writer.writerow(data)

if __name__ == "__main__":
    main()
