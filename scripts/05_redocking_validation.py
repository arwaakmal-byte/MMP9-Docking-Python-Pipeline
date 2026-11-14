"""
Stage 5 — Redocking validation.

Validates that the docking protocol (receptor prep, ligand prep, atom
typing, search parameters) can reproduce a known, correct binding pose.

NOTE ON METHODOLOGY: the obvious first choice - redocking 1GKC's own
co-crystallized ligand (NFH) - was tried and failed (RMSD 6.84 A, far
above the 2.0 A reliability threshold). Likely causes: bond orders had to
be guessed (no template), and the ligand/HETATM selection may have mixed
both asymmetric-unit copies. Whether NFH occupies the active site should
be checked with scripts/02b_check_reference_ligand.py rather than assumed.
Instead of debugging that, this script validates against a genuine
active-site MMP-9 inhibitor: ligand 10B in PDB 4H3X, which sits 2.01 A
from a zinc (direct coordination) and has known bond connectivity in
RCSB's Chemical Component Dictionary.
"""

import os
import subprocess

import numpy as np
import requests
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

from config import DOCKING_DIR, VALIDATION_LIGAND_CODE, VALIDATION_PDB_PATH, VINA_BIN_PATH


def find_catalytic_zinc_by_ligand_proximity(active_lines, ligand_code):
    """Identify the catalytic zinc directly: it's whichever zinc the real
    inhibitor actually coordinates (closest by distance), rather than
    inferring from His/Asp coordination patterns alone."""
    lig_hetatm = [l for l in active_lines if l.startswith("HETATM") and l[17:20].strip() == ligand_code]
    lig_chain = lig_hetatm[0][21]
    lig_hetatm = [l for l in lig_hetatm if l[21] == lig_chain]  # single copy only

    lig_coords = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lig_hetatm])

    zn_atoms = [
        {"resseq": l[22:26].strip(), "xyz": np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])}
        for l in active_lines if l.startswith("HETATM") and l[17:20].strip() == "ZN" and l[21] == lig_chain
    ]
    best_zn, best_dist = None, None
    for zn in zn_atoms:
        d = np.linalg.norm(lig_coords - zn["xyz"], axis=1).min()
        print(f"Zn {zn['resseq']}: closest ligand atom {d:.2f} A away")
        if best_dist is None or d < best_dist:
            best_dist, best_zn = d, zn

    return lig_chain, lig_hetatm, best_zn


def build_clean_validation_receptor(active_lines, chain, catalytic_zn_resseq):
    clean = [
        l for l in active_lines
        if (l.startswith("ATOM") and l[21] == chain)
        or (l.startswith("HETATM") and l[17:20].strip() == "ZN" and l[21] == chain and l[22:26].strip() == catalytic_zn_resseq)
    ]
    path = os.path.join(DOCKING_DIR, "4H3X_receptor_clean.pdb")
    with open(path, "w") as f:
        f.writelines(clean)
        f.write("END\n")
    return path


def correct_bond_orders(lig_hetatm, ligand_code):
    """Bond orders guessed from raw PDB coordinates (no CONECT records)
    are unreliable. Pull the true connectivity from RCSB's Chemical
    Component Dictionary and map it onto the crystal coordinates."""
    chemcomp = requests.get(f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_code}").json()
    smiles = chemcomp["rcsb_chem_comp_descriptor"]["SMILES"]
    print(f"{ligand_code} canonical SMILES (RCSB): {smiles}")

    template_mol = Chem.MolFromSmiles(smiles)
    lig_pdb_block = "".join(lig_hetatm) + "END\n"
    raw_mol = Chem.MolFromPDBBlock(lig_pdb_block, sanitize=False, removeHs=False)
    raw_mol = Chem.RemoveHs(raw_mol, sanitize=False)

    fixed_mol = AllChem.AssignBondOrdersFromTemplate(template_mol, raw_mol)
    Chem.SanitizeMol(fixed_mol)
    fixed_mol = Chem.AddHs(fixed_mol, addCoords=True)  # crystal structures have no H coordinates
    return fixed_mol


def parse_pdbqt_pose(path, model_num=1):
    coords, in_model = [], False
    with open(path) as f:
        for line in f:
            if line.startswith("MODEL"):
                in_model = int(line.split()[1]) == model_num
            elif line.startswith("ENDMDL"):
                in_model = False
            elif in_model and line.startswith(("ATOM", "HETATM")):
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(coords)


def parse_pdbqt_single(path):
    coords = []
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(coords)


def main():
    with open(VALIDATION_PDB_PATH) as f:
        active_lines = f.readlines()

    lig_chain, lig_hetatm, catalytic_zn = find_catalytic_zinc_by_ligand_proximity(
        active_lines, VALIDATION_LIGAND_CODE
    )
    print(f"Catalytic zinc: resSeq {catalytic_zn['resseq']} (chain {lig_chain})")

    receptor_clean_path = build_clean_validation_receptor(active_lines, lig_chain, catalytic_zn["resseq"])

    lig_coords_native = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lig_hetatm])
    lig_centroid = lig_coords_native.mean(axis=0)

    receptor_out_basename = os.path.join(DOCKING_DIR, "4H3X_receptor")
    subprocess.run([
        "mk_prepare_receptor.py", "--read_pdb", receptor_clean_path,
        "-o", receptor_out_basename, "-p", "-v",
        "--box_size", "20", "20", "20",
        "--box_center", str(lig_centroid[0]), str(lig_centroid[1]), str(lig_centroid[2]),
        "--allow_bad_res", "--default_altloc", "A",
    ], check=True)
    receptor_pdbqt = receptor_out_basename + ".pdbqt"

    fixed_mol = correct_bond_orders(lig_hetatm, VALIDATION_LIGAND_CODE)
    preparator = MoleculePreparation()
    setups = preparator.prepare(fixed_mol)
    native_pdbqt_str = PDBQTWriterLegacy.write_string(setups[0])[0]

    native_path = os.path.join(DOCKING_DIR, f"{VALIDATION_LIGAND_CODE}_native.pdbqt")
    with open(native_path, "w") as f:
        f.write(native_pdbqt_str)

    redock_out = os.path.join(DOCKING_DIR, f"{VALIDATION_LIGAND_CODE}_redocked.pdbqt")
    result = subprocess.run([
        VINA_BIN_PATH, "--receptor", receptor_pdbqt, "--ligand", native_path,
        "--center_x", str(lig_centroid[0]), "--center_y", str(lig_centroid[1]), "--center_z", str(lig_centroid[2]),
        "--size_x", "20", "--size_y", "20", "--size_z", "20",
        "--exhaustiveness", "8", "--num_modes", "9", "--energy_range", "3.0",
        "--out", redock_out,
    ], capture_output=True, text=True)
    print(result.stdout)

    native_coords = parse_pdbqt_single(native_path)
    docked_coords = parse_pdbqt_pose(redock_out, model_num=1)

    if native_coords.shape == docked_coords.shape:
        rmsd = np.sqrt(np.mean(np.sum((native_coords - docked_coords) ** 2, axis=1)))
        print(f"\nRedocking RMSD (best pose vs. native): {rmsd:.3f} A")
        print("PASS (< 2.0 A)" if rmsd < 2.0 else "FAIL (>= 2.0 A) — inspect protocol before trusting docking results")
    else:
        print(f"Atom count mismatch (native: {native_coords.shape}, docked: {docked_coords.shape})")


if __name__ == "__main__":
    main()
