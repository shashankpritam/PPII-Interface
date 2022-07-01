# Change the input method on the pp2_predictor.py first before running the script where input_pdb = sys.argv[1]
# If database search is required please run pp2_pred_for_database.py and provide the query PDB ID in input_database.txt
# Let the input command be as sys.argv[1]. Also please provide the PDB files in database_folder for database look-up.
import os
with open("input_database.txt") as input_pdbs:
    lines = input_pdbs.readlines()[1:]
    for pdb in lines:
        run_test = "python3 pp2_predictor.py "+pdb
        #print(run_test)
        os.system(run_test)
