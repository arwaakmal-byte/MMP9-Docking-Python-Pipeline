"""
Stage 4 — Dock all 20 ligands against the prepared receptor.

Vina's algorithm, briefly: within the grid box, it runs an iterated local
search — random starting poses refined via BFGS local optimization,
perturbed by random mutation, accepted/rejected via a Metropolis criterion
— repeated `exhaustiveness` times to reduce the risk of settling on a bad
local minimum. It returns up to `num_modes` distinct poses (clustered at
least 1.0 A RMSD apart), all within `energy_range` kcal/mol of the best
score found.
"""

import os
import subprocess

import pandas as pd

from config import (
    DOCKING_DIR,
    ENERGY_RANGE,
    EXHAUSTIVENESS,
    GRID_SIZE,
    LEAD_CHEMBL_IDS,
    NUM_MODES,
    PROJECT_DIR,
    VINA_BIN_PATH,
)


def dock_ligand(receptor_pdbqt_path, ligand_path, grid_center, out_path):
    cmd = [
        VINA_BIN_PATH,
        "--receptor", receptor_pdbqt_path,
        "--ligand", ligand_path,
        "--center_x", str(grid_center[0]),
        "--center_y", str(grid_center[1]),
        "--center_z", str(grid_center[2]),
        "--size_x", str(GRID_SIZE), "--size_y", str(GRID_SIZE), "--size_z", str(GRID_SIZE),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", str(NUM_MODES),
        "--energy_range", str(ENERGY_RANGE),
        "--out", out_path,
    ]
    return subprocess.run(cmd, capture_output=True, text=True)


def parse_best_affinity(vina_stdout):
    for line in vina_stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "1":
            try:
                return float(parts[1])
            except ValueError:
                pass
    return None


def main(grid_center):
    receptor_pdbqt_path = os.path.join(DOCKING_DIR, "receptor.pdbqt")
    results = []

    for cid in LEAD_CHEMBL_IDS:
        ligand_path = os.path.join(DOCKING_DIR, f"{cid}.pdbqt")
        out_path = os.path.join(DOCKING_DIR, f"{cid}_docked.pdbqt")
        log_path = os.path.join(DOCKING_DIR, f"{cid}_log.txt")

        print(f"Docking {cid}...")
        result = dock_ligand(receptor_pdbqt_path, ligand_path, grid_center, out_path)
        with open(log_path, "w") as f:
            f.write(result.stdout)

        affinity = parse_best_affinity(result.stdout)
        if affinity is not None:
            print(f"  Best affinity: {affinity} kcal/mol")
        else:
            print(f"  Could not parse affinity — see {log_path}")

        results.append({"compound_id": cid, "best_affinity_kcal_mol": affinity})

    df = pd.DataFrame(results).sort_values("best_affinity_kcal_mol")
    out_csv = os.path.join(PROJECT_DIR, "MMP9_docking_results.csv")
    df.to_csv(out_csv, index=False)
    print(f"\nSaved: {out_csv}")
    return df


if __name__ == "__main__":
    # Grid center = catalytic zinc coordinates. Re-read them directly from
    # the cleaned receptor PDB (written by 02_receptor_prep.py) rather than
    # hardcoding, so this stays correct if the receptor is ever rebuilt.
    import numpy as np
    from config import CATALYTIC_ZN_CHAIN, CATALYTIC_ZN_RESSEQ, RECEPTOR_PDB_ID

    receptor_clean_path = os.path.join(DOCKING_DIR, f"{RECEPTOR_PDB_ID}_chainA_receptor.pdb")
    with open(receptor_clean_path) as f:
        for line in f:
            if (line.startswith("HETATM") and line[17:20].strip() == "ZN"
                    and line[21] == CATALYTIC_ZN_CHAIN and line[22:26].strip() == CATALYTIC_ZN_RESSEQ):
                center = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                break
        else:
            raise RuntimeError("Catalytic zinc not found in cleaned receptor PDB — run 02_receptor_prep.py first.")

    main(center)
