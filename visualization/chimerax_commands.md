# ChimeraX visualization commands

Install ChimeraX (free for academic/non-commercial use): https://www.cgl.ucsf.edu/chimerax/download.html

## Quick start (recommended)

1. Run `scripts/08_make_chimerax_pack.py` (last notebook cell) and unzip `chimerax_pack.zip` into a folder.
2. In ChimeraX, run (use your real path):
   ```
   cd ~/Downloads/chimerax_pack
   open CHEMBL75164.cxc
   ```
   Each `<compound>.cxc` opens the receptor + best pose, shows the pocket (5 Å), the zinc and its histidines, a zinc-distance line, H-bonds, and saves `<compound>_pose.png`.
3. Other scripts: `show_all_top5.cxc` (overlay of the top 5), `validation_4H3X.cxc` (crystal vs redocked pose + RMSD), `grid_box.cxc` (docking box).

## Manual commands (type into the ChimeraX command line)

| Goal | Command |
|---|---|
| Open receptor and pose | `open receptor.pdbqt` then `open CHEMBL75164_best.pdbqt` |
| Open all 9 poses | `open CHEMBL75164_docked.pdbqt` (opens as sub-models), then `viewdockx` |
| Cartoon only, transparent | `hide atoms` · `cartoon` · `transparency 60 target c` |
| Ligand as sticks | `show #2 atoms` · `style #2 stick` · `color #2 byhetero` |
| Zinc as purple sphere | `show :ZN atoms` · `style :ZN sphere` · `color :ZN purple target a` |
| Pocket residues | `show (#1 & protein) & (#2 :<5) atoms` |
| Hydrogen bonds | `hbonds #2 restrict cross reveal true` |
| Atom contacts / clashes | `contacts #2 restrict cross reveal true` |
| Measure distance | `distance #1:ZN@ZN #2@O1` (use the atom name from the CSV) |
| Remove distances | `~distance` |
| Label residues | `label (#1 & protein) & (#2 :<5) residues` |
| Protein surface | `surface #1` · `transparency 70 target s` |
| RMSD of two ligands | `rmsd #3 to #2 matchAtomNames true` |
| White background | `set bgColor white` |
| Save image | `save pose.png width 1600 height 1200 supersample 3` |
| Save session | `save mmp9.cxs` |
| Help for any command | `usage hbonds` |

## Tips

- Model numbers (`#1`, `#2`) depend on open order. If a command says "no atoms", check the model panel (Tools → General → Model Panel).
- A multi-pose file opens as `#2.1 … #2.9`; use `#2.1` for the best pose, or the `*_best.pdbqt` files which hold only pose 1.
- Vina PDBQT files carry no bond-order information; ChimeraX infers bonds from geometry, so double bonds may look single. This is a display matter only.
- Vina scores are estimates (kcal/mol), not measured affinities. The zinc distance and a sensible pose matter more than small score differences.

## Command-line ChimeraX (no GUI clicking)

```
chimerax --nogui --offscreen --script "open CHEMBL75164.cxc; exit"
```
