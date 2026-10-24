# 6LU7 · SARS-CoV-2 main protease

Provenance is summarised in [../README.md](../README.md). This file explains
exactly how the raw structure becomes Fairouz's docking receptor.

## Reproduce the preparation

```bash
python -m backend.receptor
```

This regenerates everything in `prepared/` from `raw/6LU7.pdb` and rewrites
`prepared/preparation_log.json`. The raw file is checked against its SHA-256
checksum first and is never modified. Prepared files are written with LF
line endings, the way Git stores them, so the checksum recorded in the log
matches a fresh checkout on any platform.

## What the raw file contains

| Chain | Contents | Used? |
| ----- | -------- | ----- |
| A | Main protease, residues 1–306 (2,367 atoms) | Kept as the receptor |
| A (HETATM) | 84 crystallographic waters | Removed |
| C | Inhibitor N3: 02J–ALA–VAL–LEU–PJE–010 (49 atoms), covalently bonded to Cys145 SG | Removed from the receptor; used only to define the box |

There are no alternate locations, and the file has no hydrogens.

## Transformations (in order)

1. Verify the raw file's SHA-256 checksum.
2. Keep model 1 and chain A only.
3. Remove N3 (chain C). In the crystal, N3 is covalently bonded to Cys145;
   after removal, Cys145 is treated as a normal free cysteine.
4. Remove all 84 waters. Some active-site waters can matter for binding,
   but keeping a few by hand would be an arbitrary choice for a demo.
5. Remove any hydrogens, then write `6LU7_receptor_clean.pdb`.
6. Run Meeko `mk_prepare_receptor` (Meeko 0.8). It matches each residue to a
   template, adds polar hydrogens, assigns Gasteiger charges and AutoDock
   atom types, and writes the rigid receptor `6LU7_receptor.pdbqt`.
   Histidines get Meeko's default HIE state (hydrogen on Nε). This is a
   simplification: His41, in the catalytic dyad, could be protonated
   differently in reality.

## Binding box

The box is derived from N3's coordinates by
`backend.receptor.calculate_binding_box()`:

- **Centre:** the midpoint of N3's extent (minimum to maximum) on each axis.
- **Size:** N3's extent on each axis plus **5 Å of padding on each side**.
- **Units:** ångströms, in the coordinate frame of the PDB file.

| | x | y | z |
| --- | --- | --- | --- |
| N3 minimum | −14.397 | 1.458 | 63.232 |
| N3 maximum | −5.138 | 21.414 | 74.575 |
| **Centre** | **−9.7675** | **11.436** | **68.9035** |
| **Size** | **19.259** | **29.956** | **21.343** |

These values are stored in `config.py` (`BINDING_BOX_CENTER`,
`BINDING_BOX_SIZE`). A unit test recomputes them from the raw file, so the
config can't drift from the documented method.

**Why 5 Å of padding:** it gives ligands room to rotate and sample the
pocket around where N3 sits, while keeping every edge at or under about
30 Å, the range in which Vina's search is reliable at default
exhaustiveness. The padding was fixed before any compound was docked and was
not tuned to improve any compound's score.

**Limitations:** N3 is long and peptide-like, so the box is elongated along
y and covers more than the S1 pocket alone. A box based on one ligand
assumes that compounds bind where N3 binds; other binding sites are not
explored.

## Reference ligand sanity check

Re-docking a crystal ligand and checking that it returns near its crystal
pose is the standard sanity check for a docking setup. **It has not been
done for N3, deliberately**:

- N3 is a *covalent* inhibitor, and Vina docks non-covalently. A
  non-covalent re-dock of a covalent ligand is not a fair test.
- N3 is large and very flexible (6 residues), well beyond the range where
  default Vina settings reproduce poses reliably.
- Rebuilding N3's correct bond orders and hydrogens from the PDB
  coordinates is a separate task.

What *is* checked now (`python -m pytest -m integration`): Vina really runs;
returned scores are finite and negative; docked poses lie inside the box;
and the same seed gives the same score.

**Planned validation task:** re-dock a *non-covalent* co-crystallised Mpro
inhibitor (for example X77 from PDB 6W63) and measure the pose RMSD to its
crystal position. This is a software and docking-setup check, not
biological validation.
