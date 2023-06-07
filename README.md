
# PPII Binding Site Predictor

This program analyzes a PDB structure to identify possible binding sites of the PPII helix along with the best PPII template.

## Usage

This program takes a PDB structure as input and returns the possible binding site of the PPII helix.

**Input:** 4-letter PDB ID (e.g., 1CKA)

**Output:**
- If the given structure does not contain a PPII binding site, you will receive a notification.
- If the given structure does contain a PPII binding site, you will receive the structure with PPII bound to it in the "best possible orientation" at the predicted binding site of the PPII.

## Author

- Shashank Pritam (shashankpritam[at]gmail[dot]com)

## License

This project is licensed under the LGPL-2.1 License.

## Requirements

- numpy
- scipy
- modeller
- biopython
- Python 3.8.10
- Tested on WSL Ubuntu 20.4

**Note:** An active internet connection is required if PDB files are not provided in the `pp2pred` folder (database_folder).

## Usage

1. Make sure you have the required modules installed (see Requirements section).
2. Ensure you have Python version 3.8.10 installed.
3. Clone the repository or download the source code.
4. Provide an input file `input_database.txt` with a list of 4-letter PDB IDs (e.g., `1CKA`) to analyze. One PDB ID per line.
5. Execute the program by running the following command in the terminal:

`python3 pp2_pred_for_database.py input_database.txt`

6. The program will analyze each PDB ID in the input file and generate the following output:
   - If no binding site for PPII is found, a notification will be provided in the log file `pp2_pred_db_log.log`.
   - If a binding site is found, the program will generate a "CLICK" transferred PPII on the binding site, saved as `input_pdb_best_template_model_pdb_4sim.pdb`.
   - Additionally, `input_pdb_best_template_model_pdb_result.pdb` will contain all the accepted snapshots of the Monte Carlo moves.
   - The program will generate two CSV files as output: `output.csv` and `binding_site_output.csv`. These files contain the predicted binding site details for each PDB ID.
     The columns in the CSV files include: PDB_ID, Model, Chain, TRP, NBR, Template_PPII, and Score.

For any help, please contact through the provided email ID.
