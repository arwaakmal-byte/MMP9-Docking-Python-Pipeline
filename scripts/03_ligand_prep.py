"""
Stage 3 — Ligand preparation for the 20 candidate compounds.

Structures are pulled live from the ChEMBL REST API by compound ID rather
than hand-typed from memory, so a wrong or mistyped structure fails loudly
(a bad API call) instead of silently docking the wrong molecule.
"""

import os

import requests
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

from config import DOCKING_DIR, LEAD_CHEMBL_IDS


def fetch_smiles(chembl_id):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
    r = requests.get(url)
    r.raise_for_status()
    return r.json()["molecule_structures"]["canonical_smiles"]


def prepare_ligand(chembl_id):
    smi = fetch_smiles(chembl_id)
    print(f"{chembl_id}: {smi}")

    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)

    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)
    pdbqt_string = PDBQTWriterLegacy.write_string(mol_setups[0])[0]

    out_path = os.path.join(DOCKING_DIR, f"{chembl_id}.pdbqt")
    with open(out_path, "w") as f:
        f.write(pdbqt_string)
    return out_path


if __name__ == "__main__":
    prepared, failed = {}, []
    for cid in LEAD_CHEMBL_IDS:
        try:
            prepared[cid] = prepare_ligand(cid)
            print(f"  Saved {prepared[cid]}\n")
        except Exception as e:
            print(f"  Failed on {cid}: {e}\n")
            failed.append(cid)

    print(f"Prepared {len(prepared)}/{len(LEAD_CHEMBL_IDS)} ligands.")
    if failed:
        print(f"Failed: {failed}")
