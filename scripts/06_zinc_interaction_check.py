"""
Stage 6 — Zinc-coordination check for all 20 docked poses.

A strong Vina affinity is not, on its own, evidence that a compound
engages MMP-9's catalytic mechanism — it can just as easily come from
nonspecific hydrophobic packing. This script checks, for each compound's
best pose, whether a real heteroatom (O/N/S — capable of donating a lone
pair to Zn2+) sits close to the catalytic zinc, the way a genuine MMP
inhibitor's zinc-binding group (hydroxamate, sulfonamide, etc.) should.

Only heavy atoms are considered — a polar hydrogen being closest to the
zinc isn't chemically meaningful, since coordination happens through a
lone pair on the parent heteroatom, not through the hydrogen itself.
"""

import os

import numpy as np
import pandas as pd

from config import CATALYTIC_ZN_CHAIN, CATALYTIC_ZN_RESSEQ, DOCKING_DIR, LEAD_CHEMBL_IDS, PROJECT_DIR, RECEPTOR_PDB_ID

COORDINATING_TYPES = {"OA", "N", "NA", "SA", "O", "S"}
ZN_COORD_CUTOFF = 2.6  # Angstroms
ZN_PROXIMITY_CUTOFF = 4.0


def get_catalytic_zn_coord():
    receptor_clean_path = os.path.join(DOCKING_DIR, f"{RECEPTOR_PDB_ID}_chainA_receptor.pdb")
    with open(receptor_clean_path) as f:
        for line in f:
            if (line.startswith("HETATM") and line[17:20].strip() == "ZN"
                    and line[21] == CATALYTIC_ZN_CHAIN and line[22:26].strip() == CATALYTIC_ZN_RESSEQ):
                return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    raise RuntimeError("Catalytic zinc not found — run 02_receptor_prep.py first.")


def parse_pdbqt_best_pose_atoms(path):
    atoms = []
    in_model = False
    with open(path) as f:
        for line in f:
            if line.startswith("MODEL"):
                in_model = int(line.split()[1]) == 1
            elif line.startswith("ENDMDL"):
                in_model = False
            elif in_model and line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip()
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                ad_type = line[77:79].strip()
                atoms.append((atom_name, ad_type, xyz))
    return atoms


def main():
    zn_coord = get_catalytic_zn_coord()
    print(f"Catalytic zinc coordinates: {zn_coord}")

    results = []
    for cid in LEAD_CHEMBL_IDS:
        docked_path = os.path.join(DOCKING_DIR, f"{cid}_docked.pdbqt")
        atoms = parse_pdbqt_best_pose_atoms(docked_path)
        heavy_atoms = [(n, t, xyz) for n, t, xyz in atoms if t not in ("H", "HD")]

        best_dist, best_atom = None, None
        for name, ad_type, xyz in heavy_atoms:
            d = np.linalg.norm(xyz - zn_coord)
            if best_dist is None or d < best_dist:
                best_dist, best_atom = d, (name, ad_type)

        coordinating = best_atom[1] in COORDINATING_TYPES and best_dist <= ZN_COORD_CUTOFF
        results.append({
            "compound_id": cid,
            "closest_heavy_atom": best_atom[0],
            "closest_heavy_atom_type": best_atom[1],
            "dist_to_catalytic_zn": round(best_dist, 2),
            "zinc_coordinating": coordinating,
            "near_zinc_pocket": best_dist <= ZN_PROXIMITY_CUTOFF,
        })

    df = pd.DataFrame(results).sort_values("dist_to_catalytic_zn")
    out_path = os.path.join(PROJECT_DIR, "MMP9_zinc_interaction_check.csv")
    df.to_csv(out_path, index=False)
    print(f"\nSaved: {out_path}")
    print(df.to_string(index=False))
    return df


if __name__ == "__main__":
    main()
