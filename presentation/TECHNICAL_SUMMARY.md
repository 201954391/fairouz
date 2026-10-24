# Fairouz: technical summary

A one-page overview for judges. The full detail is in the repository
README.

## What it is

A computational virtual-screening research prototype. It docks compounds
into a viral protein, ranks them with a transparent score, explains every
score, and writes a reproducible report. It prioritises compounds for further
investigation. It does not establish binding, efficacy, safety or clinical
usefulness.

## Architecture

```
Streamlit UI (Home · Screening · Results)
      │
Screening pipeline: status, timing and detail for every stage
  1 Compound loading   → 2 Molecular analysis (RDKit) → 3 Ligand preparation (RDKit + Meeko)
  → 4 Docking (AutoDock Vina) → 5 Ranking → 6 Explainability → 7 Report (MD · JSON · HTML)
```

Each stage consumes the previous stage's real output; nothing is
recalculated downstream. The backend has no UI code, and the UI has no
scientific logic.

## Real scientific components

| Component | Implementation |
| --------- | -------------- |
| Target | SARS-CoV-2 main protease, PDB 6LU7 (X-ray, 2.16 Å). Raw file checksum-verified. Chain A kept; the inhibitor N3 and 84 waters removed; Meeko adds polar H and charges. |
| Binding box | Derived from N3's coordinates plus 5 Å padding. Stored in config and recomputed by a test. |
| Compounds | 10 public molecules; every SMILES matches PubChem. |
| Descriptors | RDKit: MW, Wildman–Crippen LogP, TPSA, H-bond donors and acceptors. |
| Ligands | RDKit ETKDGv3 conformer (seed 42), MMFF94, Meeko PDBQT. |
| Docking | AutoDock Vina 1.2.7 executable; exhaustiveness 8, 9 poses, seed 42. |

## Ranking approach

**Screening Priority Score = docking points (≤ 70) + property points (≤ 30).**

- **Docking:** a fixed scale (−4 kcal/mol → 0, −10 kcal/mol → 70). It
  doesn't depend on which other compounds are in the batch.
- **Properties:** the mean of five sub-scores using documented
  Lipinski/Ghose/Veber-style ranges.
- **Ties:** broken by more negative affinity, then compound ID.

There is no supervised ML: 10 molecules with no activity labels can't train
a meaningful model, so the absence of fake ML is deliberate.

## Explainability approach

Fixed rules applied to the ranking's own numbers; no language model, and no
SHAP because there is no model. For each compound it shows:

- docking and property points, which add up exactly to the score
- each descriptor's value, configured range, position and sub-score
- the main positive and limiting factors
- the primary driver
- comparisons with neighbouring compounds, and tie-break notes

Unranked compounds show the failed stage and the recorded error.

## Reproducibility

- **Fixed seeds:** two real runs give identical docking scores, ranking,
  explanations and report content. This is checked by an integration test.
- **Run ID:** UTC time + target + configuration hash + results hash.
  Different outcomes always get different IDs.
- **Reports record:** software versions, the Git commit, and SHA-256
  checksums of the dataset, raw target and prepared receptor.
- **Cross-platform files:** the prepared receptor is written with LF line
  endings, so its checksum is stable on every platform; regenerating it
  reproduces the committed file byte for byte.

## Testing

- **235 unit tests** (about 15 s): every backend module, display
  formatting, degraded modes, report serialisation, and wording scans for
  unsupported claims.
- **5 integration tests:** real Vina docking, a reproducible full pipeline
  run, and receptor regeneration.

## Responsible-use design

- Public data only; no patient data or personal data.
- Provenance is documented for every compound and for the target.
- One concise disclaimer on each result surface, with full limitations in
  the reports.
- Neutral, colour-independent status displays, and cautious wording
  ("predicted docking affinity", "screening priority").
- Failures are shown honestly: a compound without a valid docking result
  never receives a score.

## Limitations

- One rigid receptor structure, with no waters and a fixed protonation state.
- A small demo set of well-known, weakly binding molecules.
- Docking affinities are computational estimates, not measured binding.
- The ranking weights and ranges are heuristic and not validated against
  experimental data.
- The docking setup has not yet been validated by redocking a non-covalent
  reference ligand.
- Every result needs experimental validation.
