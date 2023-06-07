# Security Policy

## Supported Versions

Use this section to tell people about which versions of your project are currently being supported with security updates.

| Version | Supported          |
| ------- | ------------------ |
| 1.0.x   | :white_check_mark: |


## Reporting a Vulnerability

Use this section to tell people how to report a vulnerability.

If you discover a vulnerability in this program, please report it to us promptly. We appreciate your efforts in disclosing the issue responsibly and will work to address it as quickly as possible.

To report a vulnerability, please follow these steps:

1. Contact me at my mail ID provided below to share the details of the vulnerability.
1. Provide a clear and concise description of the vulnerability, including any relevant information or proof-of-concept code.
1. Include your contact information (name and email address) so that we can communicate with you regarding the vulnerability.
1. We will investigate the reported vulnerability and determine its impact and severity.
1. If the vulnerability is accepted, we will take appropriate measures to address it and develop a patch or mitigation strategy.
1. Once the vulnerability is resolved, we will notify you and provide any necessary updates on the fix.
1. If the vulnerability is deemed out of scope or does not pose a significant risk, we will inform you of our decision.

We are committed to ensuring the security of our program and appreciate your assistance in making it better. Thank you for helping us maintain a safe and secure environment.

---

# Program Information

This program analyzes a PDB structure to identify possible binding sites of the PPII helix along with the best PPII template.

**Input:** 4-letter PDB ID (e.g., 1CKA)

**Output:** If no binding site for PPII is found, a notification will be provided in the log file (`pp2_pred_db_log.log`). If a binding site is found, the program will generate a "CLICK" transferred PPII on the binding site, saved as `input_pdb_best_template_model_pdb_4sim.pdb`. Additionally, `input_pdb_best_template_model_pdb_result.pdb` will contain all the accepted snapshots of the Monte Carlo moves.

**Author:** Shashank Pritam ([shashankpritam@gmail.com](mailto:shashankpritam@gmail.com))

**License:** LGPL

---

## System Requirements

- Python Version: 3.8.10
- Tested on: WSL Ubuntu 20.4

---

## Dependencies

The following Python modules are required to run this program:

- numpy
- scipy
- modeller
- biopython

---

## Internet Connection

An active internet connection is required if PDB files are not provided in the `pp2pred` folder (database_folder).

