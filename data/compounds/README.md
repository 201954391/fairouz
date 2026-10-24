# Compound datasets

## What `compounds.csv` contains

Ten well-known, non-sensitive small molecules (simple chemicals and common
over-the-counter medicines), one per row:

| Column | Meaning |
| ------ | ------- |
| `compound_id` | Fairouz identifier (`FZ-0001` …), unique per row |
| `name` | Common name of the molecule |
| `smiles` | Structure as a SMILES string |
| `source` | Where the structure came from (PubChem compound ID) |

## Why it exists

The dataset exists **only to test and demonstrate the cheminformatics
pipeline**: loading, SMILES validation and RDKit descriptor calculation. The
molecules were chosen because their structures are simple, widely documented
and easy to check.

> **Inclusion in this dataset does not imply antiviral activity or any
> therapeutic effect.** None of these molecules is presented as a treatment
> for any virus.

## Provenance

Every SMILES string was retrieved from PubChem (National Center for
Biotechnology Information) on 2026-10-24 using the PUG REST service and
matches PubChem's record exactly. As an independent check, the molecular
weight RDKit calculates for every molecule matches the weight PubChem reports.

| ID | Name | PubChem CID | Note |
| -- | ---- | ----------- | ---- |
| FZ-0001 | Ethanol | 702 | |
| FZ-0002 | Benzene | 241 | |
| FZ-0003 | Urea | 1176 | |
| FZ-0004 | Glycine | 750 | |
| FZ-0005 | Benzoic acid | 243 | |
| FZ-0006 | Salicylic acid | 338 | |
| FZ-0007 | Aspirin | 2244 | |
| FZ-0008 | Paracetamol | 1983 | PubChem title: Acetaminophen |
| FZ-0009 | Ibuprofen | 3672 | Racemic, no stereochemistry |
| FZ-0010 | Caffeine | 2519 | |

To reproduce the retrieval:

```
https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/702,241,1176,750,243,338,2244,1983,3672,2519/property/Title,SMILES,MolecularWeight/CSV
```

## Licensing and source considerations

PubChem is a free, public database run by the U.S. National Library of
Medicine. Basic facts such as a molecule's name and structure are not
generally subject to copyright. PubChem notes that some depositor-supplied
data can carry its own terms, so before adding larger datasets, check the
terms for each source and record them here.

## Limitations

- Ten molecules are far too few to represent the chemical space relevant to
  drug discovery. The set is for software testing, not scientific conclusions.
- The molecules were picked by hand for simplicity, which is itself a bias.
- Descriptor definitions differ between tools: PubChem's TPSA and
  hydrogen-bond acceptor counts can differ slightly from RDKit's for the
  same molecule.
- Stereochemistry is not represented for these molecules.

## Rules for adding a dataset

- Only use publicly available data whose licence allows this use.
- Record the source, version or retrieval date, and licence in this README.
- Note known gaps, biases or quality issues.
- Never add patient or personal data.
- Keep the required columns: `compound_id`, `name`, `smiles`, `source`.
  Rows with missing values are skipped and reported by the loader, and
  duplicate IDs stop loading entirely.
