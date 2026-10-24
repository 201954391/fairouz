# Viral protein targets

Each target lives in its own folder, with the raw downloaded structure and
the prepared docking files kept separate:

```
targets/
└── 6LU7/
    ├── README.md                      Target-specific provenance and preparation
    ├── raw/6LU7.pdb                   Original RCSB file, never modified
    └── prepared/
        ├── 6LU7_receptor_clean.pdb    Protein chain A, heavy atoms only
        ├── 6LU7_receptor.pdbqt        Rigid receptor used by AutoDock Vina
        └── preparation_log.json       Every step, checksum, box and tool version
```

## Current target: SARS-CoV-2 main protease (PDB 6LU7)

| Field | Value |
| ----- | ----- |
| Protein | 3C-like proteinase (main protease, Mpro, 3CLpro, nsp5) |
| Organism | Severe acute respiratory syndrome coronavirus 2 (taxonomy ID 2697049) |
| PDB ID | 6LU7 (DOI 10.2210/pdb6LU7/pdb) |
| Source | RCSB Protein Data Bank, https://www.rcsb.org/structure/6LU7 |
| Method | X-ray diffraction, 2.16 Å resolution (R-free 0.235) |
| Bound ligand | Peptidomimetic inhibitor N3, covalently bonded to Cys145 |
| Released | 2020-02-05 |
| Downloaded | 2026-10-24 from `https://files.rcsb.org/download/6LU7.pdb` |
| SHA-256 (raw) | `385497151c2c093d764ffbc694589c2c1db35be638e49d715d4d8af38903966e` |
| Citation | Jin, Z. *et al.* Structure of Mpro from SARS-CoV-2 and discovery of its inhibitors. *Nature* **582**, 289–293 (2020). doi:10.1038/s41586-020-2223-y |

### Why this target

- The main protease is essential for the virus's replication and is one of
  the most studied SARS-CoV-2 targets, so its binding site is well
  characterised.
- 6LU7 was one of the first SARS-CoV-2 Mpro crystal structures and is
  widely used in teaching and benchmarking.
- It contains a **known bound inhibitor (N3)** in the active site. Fairouz
  uses N3's coordinates to place the docking box, instead of choosing the
  region by hand.

### Responsible use

Fairouz uses this structure **only as a computational demonstration** of a
docking workflow. Docking a compound into 6LU7 does not show that the
compound inhibits the protease or has any antiviral effect. Every docking
result requires experimental validation.

PDB data are freely available under the wwPDB's CC0 1.0 public domain
dedication; citing the original publication is expected.

See [6LU7/README.md](6LU7/README.md) for how the structure is prepared and
how the binding box is derived.

## Rules for adding a target

- Use publicly available structures and record the PDB ID, method,
  resolution, citation, download date and checksum.
- Keep the raw file unchanged in `raw/`, and write all prepared files to
  `prepared/` from a script, never by hand.
- Document every preparation step and how the binding site was defined.
