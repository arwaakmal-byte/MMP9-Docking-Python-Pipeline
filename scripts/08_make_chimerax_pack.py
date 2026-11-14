"""
Stage 8 - Build a ChimeraX visualization pack.

For every docked compound this writes:
  * <id>_best.pdbqt  - ONLY the best pose (model 1), so ChimeraX opens a
                       single clean model instead of 9 poses
  * <id>.cxc         - ChimeraX command script: opens receptor + pose,
                       styles them, shows the pocket, draws a zinc-distance
                       line, H-bonds, and saves a PNG
plus:
  * show_all_top5.cxc, validation_4H3X.cxc, grid_box.cxc
  * chimerax_pack.zip - everything above (download this, unzip next to
    receptor.pdbqt, then in ChimeraX:  open show_all_top5.cxc)

Zinc-distance lines use the closest ligand HEAVY atom recorded in
MMP9_zinc_interaction_check.csv, so they match the numbers in the report.
"""

import os
import shutil
import zipfile

import pandas as pd

from config import DOCKING_DIR, LEAD_CHEMBL_IDS, PROJECT_DIR, VALIDATION_LIGAND_CODE, VALIDATION_PDB_ID

PACK_DIR = os.path.join(PROJECT_DIR, "chimerax_pack")


def best_pose_only(src, dst):
    """Keep the contents of MODEL 1 only (or the whole file if no MODEL records)."""
    out, in_model, seen_model = [], True, False
    with open(src) as f:
        for line in f:
            if line.startswith("MODEL"):
                seen_model = True
                in_model = line.split()[1] == "1"
            elif line.startswith("ENDMDL"):
                if in_model:
                    break
            elif in_model and line.startswith(("ATOM", "HETATM", "ROOT", "ENDROOT",
                                               "BRANCH", "ENDBRANCH", "TORSDOF")):
                out.append(line)
    with open(dst, "w") as f:
        f.writelines(out)


def compound_cxc(cid, atom_name, dist, affinity):
    return f"""# {cid}: best pose, {affinity} kcal/mol, closest heavy atom {atom_name} to Zn = {dist} A
close session
open receptor.pdbqt
open {cid}_best.pdbqt
hide #1,2 atoms
cartoon #1
color #1 lightgray target c
transparency #1 60 target c
show #2 atoms
style #2 stick
color #2 byhetero
color #2 orange target c
# catalytic zinc and the residues coordinating it
show #1:ZN atoms
style #1:ZN sphere
color #1:ZN purple target a
show (#1 & protein) & (#1:ZN :<3) atoms
# pocket residues within 5 A of the ligand
show (#1 & protein) & (#2 :<5) atoms
style (#1 & protein) & (#2 :<5) stick
# zinc distance line to the ligand atom used in the report
distance #1:ZN@ZN #2@{atom_name}
distance style color black radius 0.12 decimalPlaces 2
# polar contacts
hbonds #2 restrict cross reveal true color red radius 0.1 dashes 6
view #2
lighting soft
graphics silhouettes true
save {cid}_pose.png width 1600 height 1200 supersample 3
"""


def main():
    os.makedirs(PACK_DIR, exist_ok=True)
    zn_df = pd.read_csv(os.path.join(PROJECT_DIR, "MMP9_zinc_interaction_check.csv")).set_index("compound_id")
    dock_df = pd.read_csv(os.path.join(PROJECT_DIR, "MMP9_docking_results.csv")).set_index("compound_id")
    final_path = os.path.join(PROJECT_DIR, "MMP9_FINAL_lead_selection.csv")
    order = (pd.read_csv(final_path)["compound_id"].tolist()
             if os.path.exists(final_path) else list(LEAD_CHEMBL_IDS))

    shutil.copy(os.path.join(DOCKING_DIR, "receptor.pdbqt"), PACK_DIR)
    for cid in order:
        best_pose_only(os.path.join(DOCKING_DIR, f"{cid}_docked.pdbqt"),
                       os.path.join(PACK_DIR, f"{cid}_best.pdbqt"))
        row = zn_df.loc[cid]
        with open(os.path.join(PACK_DIR, f"{cid}.cxc"), "w") as f:
            f.write(compound_cxc(cid, row["closest_heavy_atom"], row["dist_to_catalytic_zn"],
                                 dock_df.loc[cid, "best_affinity_kcal_mol"]))

    # Overlay of the final top 5 on one receptor
    colors = ["orange", "cyan", "magenta", "yellow", "lime"]
    lines = ["close session", "open receptor.pdbqt", "cartoon #1",
             "color #1 lightgray target c", "transparency #1 70 target c",
             "show #1:ZN atoms", "style #1:ZN sphere", "color #1:ZN purple target a"]
    for i, cid in enumerate(order[:5], start=2):
        lines += [f"open {cid}_best.pdbqt", f"style #{i} stick",
                  f"color #{i} byhetero", f"color #{i} {colors[i-2]} target c"]
    lines += ["show (#1 & protein) & (#2-6 :<5) atoms", "view #2-6", "lighting soft",
              "graphics silhouettes true", "save top5_overlay.png width 1600 height 1200 supersample 3"]
    open(os.path.join(PACK_DIR, "show_all_top5.cxc"), "w").write("\n".join(lines) + "\n")

    open(os.path.join(PACK_DIR, "validation_4H3X.cxc"), "w").write(f"""# Redocking validation: crystal pose (green) vs redocked pose (orange)
close session
open {VALIDATION_PDB_ID}_receptor_clean.pdbqt
open {VALIDATION_LIGAND_CODE}_native.pdbqt
open {VALIDATION_LIGAND_CODE}_redocked_best.pdbqt
cartoon #1
color #1 lightgray target c
transparency #1 60 target c
style #2,3 stick
color #2 lime target c
color #3 orange target c
show (#1 & protein) & (#2,3 :<5) atoms
rmsd #3 to #2 matchAtomNames true
view #2,3
lighting soft
save validation_overlay.png width 1600 height 1200 supersample 3
""")
    # best pose of the redocked validation ligand, if present
    red = os.path.join(DOCKING_DIR, f"{VALIDATION_LIGAND_CODE}_redocked.pdbqt")
    if os.path.exists(red):
        best_pose_only(red, os.path.join(PACK_DIR, f"{VALIDATION_LIGAND_CODE}_redocked_best.pdbqt"))
    for name in (f"{VALIDATION_PDB_ID}_receptor_clean.pdbqt", f"{VALIDATION_LIGAND_CODE}_native.pdbqt"):
        src = os.path.join(DOCKING_DIR, name)
        if os.path.exists(src):
            shutil.copy(src, PACK_DIR)

    # Grid box (receptor.box.pdb is written by 02_receptor_prep.py)
    open(os.path.join(PACK_DIR, "grid_box.cxc"), "w").write(
        "open receptor.pdbqt\nopen receptor.box.pdb\nstyle #2 stick\ncolor #2 red\n"
        "cartoon #1\ntransparency #1 70 target c\nview\n")
    box = os.path.join(DOCKING_DIR, "receptor.box.pdb")
    if os.path.exists(box):
        shutil.copy(box, PACK_DIR)

    zip_path = os.path.join(PROJECT_DIR, "chimerax_pack.zip")
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
        for fn in sorted(os.listdir(PACK_DIR)):
            z.write(os.path.join(PACK_DIR, fn), fn)
    print(f"Wrote {len(os.listdir(PACK_DIR))} files to {PACK_DIR}\nZip: {zip_path}")
    return zip_path


if __name__ == "__main__":
    main()
