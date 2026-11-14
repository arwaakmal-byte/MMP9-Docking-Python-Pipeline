"""
Stage 2 — Catalytic zinc identification and receptor preparation.

MMP-9's catalytic domain carries two zinc ions: one catalytic (essential
for the hydrolysis mechanism) and one structural (stabilizes the fold).
Only the catalytic zinc should be used to center the docking grid box.

Rather than rely on the co-crystallized ligand's position to indicate the
active site (see 02b_check_reference_ligand.py), the catalytic zinc is identified directly
from coordination chemistry: it is coordinated by exactly three
histidines and no aspartate, distinguishing it from the structural zinc
(three histidines PLUS one aspartate).
"""

import os
import subprocess

import numpy as np
import requests

from config import (
    CATALYTIC_ZN_CHAIN,
    CATALYTIC_ZN_RESSEQ,
    DOCKING_DIR,
    GRID_SIZE,
    RECEPTOR_PDB_ID,
    RECEPTOR_PDB_PATH,
)


def download_receptor():
    if os.path.exists(RECEPTOR_PDB_PATH):
        print(f"{RECEPTOR_PDB_ID}.pdb already present.")
        return
    resp = requests.get(f"https://files.rcsb.org/download/{RECEPTOR_PDB_ID}.pdb")
    resp.raise_for_status()
    with open(RECEPTOR_PDB_PATH, "w") as f:
        f.write(resp.text)
    print(f"Downloaded {RECEPTOR_PDB_ID}.pdb")


def parse_atom_records(lines):
    atoms = []
    for l in lines:
        if l.startswith("ATOM"):
            atoms.append({
                "chain": l[21], "resseq": l[22:26].strip(), "resname": l[17:20].strip(),
                "atom": l[12:16].strip(),
                "xyz": (float(l[30:38]), float(l[38:46]), float(l[46:54])),
            })
    return atoms


def parse_hetatm(lines, res_name):
    out = []
    for l in lines:
        if l.startswith("HETATM") and l[17:20].strip() == res_name:
            out.append({
                "chain": l[21], "resseq": l[22:26].strip(),
                "xyz": (float(l[30:38]), float(l[38:46]), float(l[46:54])),
            })
    return out


def identify_catalytic_zinc(lines):
    """Print the histidine-coordination signature for every zinc found,
    so the catalytic one can be confirmed (exactly 3 unique His, no Asp)."""
    zn_atoms = parse_hetatm(lines, "ZN")
    protein_atoms = parse_atom_records(lines)

    print(f"Total zinc atoms found: {len(zn_atoms)}")
    for zn in zn_atoms:
        zn_xyz = np.array(zn["xyz"])
        seen = {}
        for a in protein_atoms:
            if a["chain"] != zn["chain"]:
                continue
            dist = np.linalg.norm(np.array(a["xyz"]) - zn_xyz)
            if dist <= 2.8:
                key = (a["resname"], a["resseq"])
                seen.setdefault(key, dist)
        print(f"\nChain {zn['chain']}, Zn resSeq {zn['resseq']}: "
              f"{sum(1 for k in seen if k[0] == 'HIS')} unique His, "
              f"{sum(1 for k in seen if k[0] == 'ASP')} unique Asp")
        for (resname, resseq), dist in seen.items():
            print(f"  {resname}{resseq} ({dist:.2f} Å)")


def build_clean_receptor(lines):
    """Keep only Chain A protein atoms + the catalytic zinc. Waters and
    the original co-crystallized ligand (NFH) are dropped; check
    02b_check_reference_ligand.py for where NFH actually sits."""
    clean_lines = []
    zn_coord = None
    for l in lines:
        if l.startswith("ATOM") and l[21] == "A":
            clean_lines.append(l)
        elif (l.startswith("HETATM") and l[17:20].strip() == "ZN"
              and l[21] == CATALYTIC_ZN_CHAIN and l[22:26].strip() == CATALYTIC_ZN_RESSEQ):
            clean_lines.append(l)
            zn_coord = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])

    receptor_clean_path = os.path.join(DOCKING_DIR, f"{RECEPTOR_PDB_ID}_chainA_receptor.pdb")
    with open(receptor_clean_path, "w") as f:
        f.writelines(clean_lines)
        f.write("END\n")

    print(f"Cleaned receptor written: {receptor_clean_path}")
    print(f"Catalytic zinc coordinates (grid center): {zn_coord}")
    return receptor_clean_path, zn_coord


def receptor_to_pdbqt(receptor_clean_path, grid_center):
    out_basename = os.path.join(DOCKING_DIR, "receptor")
    cmd = [
        "mk_prepare_receptor.py", "--read_pdb", receptor_clean_path,
        "-o", out_basename, "-p", "-v",
        "--box_size", str(GRID_SIZE), str(GRID_SIZE), str(GRID_SIZE),
        "--box_center", str(grid_center[0]), str(grid_center[1]), str(grid_center[2]),
        "--allow_bad_res", "--default_altloc", "A",
    ]
    subprocess.run(cmd, check=True)
    receptor_pdbqt_path = out_basename + ".pdbqt"
    assert os.path.exists(receptor_pdbqt_path), "Receptor PDBQT was not written."
    print(f"Receptor PDBQT: {receptor_pdbqt_path}")
    return receptor_pdbqt_path


if __name__ == "__main__":
    download_receptor()
    with open(RECEPTOR_PDB_PATH) as f:
        lines = f.readlines()

    identify_catalytic_zinc(lines)
    receptor_clean_path, zn_coord = build_clean_receptor(lines)
    receptor_to_pdbqt(receptor_clean_path, zn_coord)
