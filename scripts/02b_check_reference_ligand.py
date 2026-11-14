"""
Stage 2b - Sanity-check the co-crystallized ligand (NFH) in 1GKC.

Reports, per chain, how many NFH atoms there are and how close the nearest
one gets to the catalytic zinc. Use this to decide whether NFH is a
trustworthy active-site reference. Also prints header records (resolution,
mutations via SEQADV, HETNAM) so you can see what the structure actually is.
"""

import numpy as np

from config import CATALYTIC_ZN_CHAIN, CATALYTIC_ZN_RESSEQ, RECEPTOR_PDB_PATH


def xyz(line):
    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])


def main():
    with open(RECEPTOR_PDB_PATH) as f:
        lines = f.readlines()

    for l in lines:
        if l.startswith(("REMARK   2 RESOLUTION", "COMPND", "SEQADV", "HETNAM")):
            print(l.rstrip())

    zn = next(xyz(l) for l in lines
              if l.startswith("HETATM") and l[17:20].strip() == "ZN"
              and l[21] == CATALYTIC_ZN_CHAIN and l[22:26].strip() == CATALYTIC_ZN_RESSEQ)

    print("\nNFH per chain (never mix chains - they are separate copies):")
    for chain in sorted({l[21] for l in lines if l.startswith("HETATM") and l[17:20].strip() == "NFH"}):
        atoms = [xyz(l) for l in lines
                 if l.startswith("HETATM") and l[17:20].strip() == "NFH" and l[21] == chain]
        d = min(np.linalg.norm(a - zn) for a in atoms)
        print(f"  chain {chain}: {len(atoms)} atoms, closest to catalytic Zn "
              f"{CATALYTIC_ZN_CHAIN}{CATALYTIC_ZN_RESSEQ}: {d:.2f} A")
    print("\nInterpretation: <~3 A from the zinc => genuine active-site ligand. "
          "Far away => not a valid redocking reference.")


if __name__ == "__main__":
    main()
