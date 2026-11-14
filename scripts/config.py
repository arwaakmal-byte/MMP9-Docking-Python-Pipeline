"""
Shared paths and constants for the MMP-9 docking pipeline.

PROJECT_DIR is expected to already contain the output files from the
upstream QSAR pipeline (composite ranking CSV, ADMET results, etc.).
Adjust it to wherever your project's persistent storage lives (e.g. a
mounted Google Drive path in Colab).
"""

import os

PROJECT_DIR = os.environ.get(
    "MMP9_PROJECT_DIR", "/content/drive/MyDrive/qsar_pdb_project"
)
DOCKING_DIR = os.path.join(PROJECT_DIR, "docking")

os.makedirs(DOCKING_DIR, exist_ok=True)

# Receptor structure used for the main 20-compound docking run
RECEPTOR_PDB_ID = "1GKC"
RECEPTOR_PDB_PATH = os.path.join(DOCKING_DIR, f"{RECEPTOR_PDB_ID}.pdb")
VINA_BIN_PATH = os.path.join(DOCKING_DIR, "vina")

# Catalytic zinc in 1GKC, Chain A — identified via 3-histidine coordination
# (His401, His405, His411; no aspartate partner, distinguishing it from the
# structural zinc at residue 1451). See scripts/02_receptor_prep.py.
CATALYTIC_ZN_CHAIN = "A"
CATALYTIC_ZN_RESSEQ = "1450"

# Docking grid box
GRID_SIZE = 20  # Angstroms per side
EXHAUSTIVENESS = 8
NUM_MODES = 9
ENERGY_RANGE = 3.0

# The 20 candidate compounds, in composite-rank order (see the upstream
# QSAR pipeline's Step 17 output: MMP9_composite_ranking_full.csv)
LEAD_CHEMBL_IDS = [
    "CHEMBL75164", "CHEMBL5285908", "CHEMBL105422", "CHEMBL309329",
    "CHEMBL106016", "CHEMBL5558206", "CHEMBL5288452", "CHEMBL1801428",
    "CHEMBL2041308", "CHEMBL1801051", "CHEMBL323371", "CHEMBL399565",
    "CHEMBL1801429", "CHEMBL400285", "CHEMBL399566", "CHEMBL78145",
    "CHEMBL2425938", "CHEMBL102504", "CHEMBL2425937", "CHEMBL319040",
]

# Validation structure: a genuine active-site MMP-9 inhibitor complex,
# used to confirm the docking protocol reproduces a known correct pose
VALIDATION_PDB_ID = "4H3X"
VALIDATION_LIGAND_CODE = "10B"
VALIDATION_PDB_PATH = os.path.join(
    PROJECT_DIR, "MMP9", "active", f"{VALIDATION_PDB_ID}.pdb"
)
