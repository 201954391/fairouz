# Fairouz

**Transparent Computational Virtual Screening** · research prototype

Fairouz is a transparent virtual-screening platform that combines RDKit
molecular descriptors, AutoDock Vina docking, deterministic ranking and
explainability to prioritise compounds for further investigation.

> Fairouz is a research and educational proof-of-concept. It is not intended
> for clinical use, medical diagnosis, treatment decisions, or direct drug
> development decisions without expert scientific review and experimental
> validation.

## At a Glance

| Question | Answer |
| -------- | ------ |
| What does it do? | Docks a set of compounds into a viral protein, scores each one transparently, ranks them, explains every score and writes a reproducible report. |
| What is real? | A real X-ray structure (SARS-CoV-2 main protease, PDB 6LU7), real RDKit descriptor calculations and real AutoDock Vina docking. Nothing is simulated or invented. |
| What is it not? | Not evidence that any compound works. Outputs are computational and need experimental validation. |
| Is there AI/ML? | No trained model. Fairouz currently uses transparent deterministic prioritisation because the available compound set does not provide an appropriate labelled training dataset. |
| What data? | 10 public molecules with PubChem-verified structures, and the PDB 6LU7 structure. No patient data and no personal data. |
| How is it tested? | 235 fast unit tests, plus 5 integration tests that run the real Vina executable and regenerate the receptor. |

**Current real result (6LU7, 10 compounds, ~23 s):**

| Rank | Compound | Score | Predicted affinity |
| ---- | -------- | ----- | ------------------ |
| 1 | Ibuprofen | 54.2 | −6.1 kcal/mol |
| 2 | Paracetamol | 43.9 | −5.2 kcal/mol |
| 3 | Aspirin | 42.7 | −5.1 kcal/mol |

These are small, well-known demo molecules with weak predicted affinities.
The point is the transparent, reproducible workflow, not the molecules.

### Quick start

```bash
pip install -r requirements.txt
```

Then install AutoDock Vina (see [Installation](#installation)) and run:

```bash
streamlit run UI/app.py
```

### Architecture

```
Streamlit UI  (UI/: Home · Screening · Results)
      │
      ▼
Screening pipeline  (backend/screening.py: status, timing and detail per stage)
  ├── 1. Compound loading      backend/compounds.py     CSV checks, PubChem provenance
  ├── 2. Molecular analysis    backend/descriptors.py   RDKit descriptors
  ├── 3. Ligand preparation    backend/ligands.py       RDKit 3D + Meeko PDBQT
  ├── 4. Docking               backend/docking.py       AutoDock Vina into PDB 6LU7
  ├── 5. Ranking               backend/ranking.py       Screening Priority Score
  ├── 6. Explainability        backend/explainability.py  rule-based breakdowns
  └── 7. Report generation     backend/report.py        Markdown · JSON · HTML
Target preparation (run once): backend/receptor.py
```

Each stage consumes the previous stage's actual output; nothing is
recalculated downstream.

### Demo workflow

1. Launch Fairouz (`streamlit run UI/app.py`). Home shows the target, the
   workflow and what is real.
2. Open **Screening**, keep 10 compounds, and press **Start screening**
   (about 25 s).
3. Open **Results → Overview** to see the top results and the stage statuses.
4. Open **Explainability** and pick a compound to see why it ranked where
   it did.
5. Open **Ranking Method** to see the formula and per-compound points.
6. Open **Report** and download the Markdown, JSON or HTML report.

The live-demo checklist and backup plan are in
[presentation/DEMO_CHECKLIST.md](presentation/DEMO_CHECKLIST.md); a short
summary for judges is in
[presentation/TECHNICAL_SUMMARY.md](presentation/TECHNICAL_SUMMARY.md).

## Overview

A researcher selects a viral protein target, and Fairouz screens a library of
chemical compounds against it. It docks each compound into the protein and
ranks the compounds with a transparent Screening Priority Score, explaining
what drove each score.

The goal is **prioritisation**: helping researchers decide which compounds
may be worth investigating further in the laboratory. Fairouz does not
identify cures, does not establish clinical effectiveness and does not
replace laboratory testing.

## Problem

Testing large numbers of compounds experimentally takes significant time,
money and laboratory resources. Many tested compounds show no meaningful
activity against the target, so a large share of that effort does not lead
to useful candidates.

## Proposed Solution

Fairouz uses computation to narrow the search before laboratory work begins.
It combines molecular property analysis and molecular docking into a ranked
shortlist, together with a transparent breakdown of every score. Researchers
can then focus experimental effort on a smaller, better-justified set of
candidates. The aim is to reduce unnecessary screening effort, not to
replace experiments.

## How It Works

```
Viral Protein Target
  → Compound Library
  → Molecular Descriptors
  → 3D Ligand Preparation
  → Molecular Docking
  → Screening Priority Ranking
  → Top Candidates
  → Research Report
  → Laboratory Validation
```

| Stage | Approach | Status |
| ----- | -------- | ------ |
| Target | SARS-CoV-2 main protease, PDB 6LU7, prepared reproducibly | **Implemented** |
| Compound library | Public compound dataset with documented provenance | **Implemented** (10-molecule demo set) |
| Molecular descriptors | RDKit: molecular weight, LogP, TPSA, H-bond donors/acceptors | **Implemented** |
| 3D ligand preparation | RDKit ETKDGv3 conformer + MMFF94 optimisation, Meeko PDBQT | **Implemented** |
| Molecular docking | AutoDock Vina 1.2: predicted docking affinity (kcal/mol) | **Implemented** |
| Priority ranking | Transparent evidence-weighted Screening Priority Score (docking 70%, properties 30%) | **Implemented** (heuristic, not ML) |
| Report | Research-style summary with limitations and uncertainty | Planned |
| Laboratory validation | Carried out by researchers; outside the software | — |

### The docking workflow in detail

1. **Target preparation** (`backend/receptor.py`, run once):
   `raw/6LU7.pdb` (checksum-verified, never modified) → keep protease chain A
   → remove the N3 inhibitor and 84 waters → `6LU7_receptor_clean.pdb`
   → Meeko adds polar hydrogens, charges and atom types → `6LU7_receptor.pdbqt`.
2. **Binding box** (`backend/receptor.py`): derived from the coordinates of
   the co-crystallised inhibitor N3, with 5 Å of padding on each side.
   Stored in `config.py`; a test recomputes it from the raw file.
3. **Ligand preparation** (`backend/ligands.py`, per compound): SMILES →
   RDKit molecule → add hydrogens → 3D conformer (ETKDGv3, seed 42) → MMFF94
   optimisation (UFF fallback) → Meeko PDBQT.
4. **Docking** (`backend/docking.py`, per compound): runs the Vina
   executable with the box, exhaustiveness 8, 9 poses and seed 42, then
   reads the affinities directly from Vina's output file.
5. **Storage**: each run gets a new UTC-timestamped folder in `results/docking/`
   containing `docking_results.csv`, `run_metadata.json` (all settings,
   checksums and tool versions), and the ligand and pose files.

6. **Ranking** (`backend/ranking.py`): every successfully docked compound gets
   a Screening Priority Score (see below), and compounds are sorted from
   highest to lowest.
7. **Explainability** (`backend/explainability.py`): each compound gets a
   rule-based breakdown of how its score was calculated (see below).
8. **Report** (`backend/report.py`): the finished results are written to
   Markdown, JSON and HTML reports in `results/reports/` (see "Report
   Generation" below).

A compound that fails at any step is marked as failed with the reason. It
never stops the batch and never receives a score.

The exact Vina command for each ligand:

```
vina --receptor data/targets/6LU7/prepared/6LU7_receptor.pdbqt \
     --ligand <compound>.pdbqt --out <compound>_poses.pdbqt \
     --center_x -9.7675 --center_y 11.436 --center_z 68.9035 \
     --size_x 19.259 --size_y 29.956 --size_z 21.343 \
     --exhaustiveness 8 --num_modes 9 --seed 42
```

## Screening Priority Score

The Screening Priority Score (0–100) ranks docked compounds **for further
computational or laboratory investigation**. It is not a probability of
efficacy, of success or of clinical benefit, and it is not a clinical
prediction.

### Inputs

All inputs are real outputs of the pipeline: the Vina docking affinity,
plus RDKit's molecular weight, LogP, TPSA, H-bond donor and H-bond acceptor
counts.

### Formula

```
Screening Priority Score = 100 × (0.7 × docking component + 0.3 × property component)
```

**Docking component (70%)** is the strongest signal. Each affinity is
placed on a fixed scale: −4 kcal/mol or weaker scores 0, −10 kcal/mol or
stronger scores 1, and values in between are linear. More negative
affinities therefore always contribute more.

We use a *fixed* scale instead of rescaling within each batch (min-max).
With min-max, the best of ten weak binders would get a perfect docking
component, and a compound's score would change depending on which other
compounds were screened alongside it. The fixed scale has neither problem,
and it cannot divide by zero when all affinities are equal. The scale's end
points are a convention covering the range where Vina scores for small
molecules typically fall. They are not calibrated against experimental data.

**Property component (30%)** is the average of five sub-scores. Each
property scores 1 inside a typical screening range, falls linearly to 0
outside it, and scores 0 beyond the outer limit:

| Property | 0 at or below | Full score | 0 at or above | Basis for the full-score range |
| -------- | ------------- | ---------- | ------------- | ------------------------------ |
| Molecular weight (g/mol) | 80 | 160–500 | 700 | Ghose (≥160), Lipinski (≤500) |
| LogP | −2.0 | −0.4 to 5.0 | 7.0 | Ghose (≥−0.4), Lipinski (≤5) |
| TPSA (Å²) | — | 0–140 | 200 | Veber (≤140) |
| H-bond donors | — | 0–5 | 8 | Lipinski (≤5) |
| H-bond acceptors | — | 0–10 | 15 | Lipinski (≤10) |

These rules of thumb describe typical oral drugs. They are a
**heuristic**: a molecule inside every range is not thereby suitable or
effective, and one outside them is not thereby useless. The fall-off widths
are a documented judgement call. All values live in `config.py`.

### Explanations

Each ranked compound gets a short rule-based sentence, for example *"Weak
predicted docking affinity (−5.2 kcal/mol); low molecular weight (outside
typical screening ranges) reduces its overall priority."* Docking is
described as strong (≤ −8 kcal/mol), moderate (≤ −6) or weak. No language
model is involved. The full breakdown is in the Explainability tab (see
"Explainability" below).

### Unranked compounds

Invalid molecules, compounds that were not docked, ligand-preparation
failures, docking failures and compounds with missing data are listed with
an explicit status and **no score**.

### Ties and reproducibility

There is no randomness. Equal scores are ordered by more negative docking
affinity, then by compound ID. Identical inputs always produce an identical
ranking, which an integration test checks end to end.

### Why there is no machine-learning model yet

Fairouz currently has 10 compounds and no measured activity data for this
target. A supervised model (random forest, neural network, and so on)
trained on that would memorise noise, and calling its output an "AI
prediction" would be misleading. Step 4 therefore deliberately uses a
transparent formula.

A supervised model may be added later **only** with an appropriate labelled
dataset, for example published, licensed Mpro activity measurements with
documented provenance. It would also need a held-out evaluation and
reported uncertainty, and it would live alongside, not silently replace,
this heuristic.

## Explainability

The **Explainability** tab shows *how Fairouz calculated* each compound's
ranking. It does not explain biological efficacy, safety, toxicity or
mechanism of action.

### How contribution values are calculated

The ranking module (`backend/ranking.py`) calculates each compound's
**docking points** (docking component × 70) and **property points**
(property component × 30). The Screening Priority Score is defined as their
exact sum. The explainability layer (`backend/explainability.py`) only
*reads* these numbers and never recalculates them with a different formula.
A test checks that the displayed contributions equal the ranking values
exactly.

For example, ibuprofen's real result: score 54.2 = docking 24.2 / 70 +
properties 30.0 / 30.

Each property's share is its sub-score × 6 (30 points ÷ 5 properties).
"Points not earned" is the maximum minus the points earned.

### How descriptor interpretations work

For each of the five descriptors, the panel shows the RDKit value, the
configured full-score range from `config.PROPERTY_SCORING_RANGES`, the
sub-score, and its position: **within**, **below** or **above the
configured range**. A value exactly on a range boundary counts as within.
The wording describes position relative to a configured screening range.
It never uses "safe/unsafe" language, and colours are neutral.

### How strengths and limitations are selected

All selections are rule-based and use only the contribution numbers:

- **Main positive factor:** the larger of the two contributions. If
  property points are at least as large as docking points, it is the
  property profile ("all five descriptors within configured ranges", or how
  many are); otherwise it is the docking affinity.
- **Main limiting factor:** the single factor that left the most points
  unearned, comparing docking (out of 70) with each property (out of 6).
  Equal losses keep the formula order: docking first, then properties.
- **Strengths:** docking counts as a strength when its component is at
  least 0.5 (about −7 kcal/mol). Properties within their ranges are also
  listed.
- **Limitations:** every factor that lost points, largest loss first.
- **Primary driver:** "Docking" if docking supplied at least 60% of the
  compound's points, "Molecular properties" if it supplied at most 40%, and
  "Both" otherwise. This describes where the points came from. It is not
  causal evidence.
- **Comparisons:** each compound is compared with its neighbours above and
  below, using the gap in docking points and property points. Gaps under
  1 point count as "similar".
- **Ties:** equal scores are explained by the tie-break rule: more negative
  docking affinity, then compound ID as the final deterministic tie-break.
  Scores that match only after rounding to one decimal place are explained
  as such.
- **Unranked compounds:** these show their status, the pipeline stage that
  stopped them (SMILES validation, 3D ligand preparation, docking, or the
  ranking input check), the recorded error, and the fact that no score was
  calculated or fabricated.

All wording thresholds live in `config.py` (`EXPLANATION_*`). They change
only the wording, never a score.

### Deterministic, and no ML explainability

Explanations are built from fixed rules. There is no language model and no
randomness, so the same ranking always produces the same explanation
(tested, including end to end with real Vina). Because there is no trained
model, ML explainability tools such as SHAP do not apply; the transparent
formula is itself the explanation.

### What the numbers do not mean

Contribution bars show points earned out of each component's maximum. They
are **not** confidence intervals, probabilities or measures of certainty.

An explanation can only be as sound as the ranking it describes: it
inherits every docking limitation and every assumption in the weights,
scale and ranges.

## Report Generation

Every screening run ends with a research report built **deterministically
from the real pipeline output**. It uses fixed templates; there is no
language model, and it recalculates nothing. Descriptor values, docking
affinities, ranking points and explanations are copied from the run exactly
as the pipeline produced them (tests check this, including with the
ranking functions disabled).

### How to generate a report

Run a screening, either in the app (**Screening → Start Screening**) or from
Python:

```bash
python -c "import config; from backend.screening import run_screening; print(run_screening(config.TARGET_LABEL, 10)['report']['files'])"
```

To regenerate a report from an existing response, call
`backend.report.generate_report(response)`.

### Formats and output folder

Reports are written to `results/reports/` in three formats sharing one name,
for example `fairouz_screening_6LU7_2026-10-24_123000.md`:

| Format | Use |
| ------ | --- |
| Markdown (`.md`) | Readable report with headings and tables (also shown in the app) |
| JSON (`.json`) | Machine-readable data for comparing runs (schema version `1.0`) |
| HTML (`.html`) | Standalone page for sharing or printing; no external resources |

- **Filenames:** Windows-safe, timestamped in UTC.
- **No overwriting:** an existing report is never overwritten; a `_2`, `_3`,
  … suffix is added instead.
- **Not committed:** `results/` is ignored by Git, so generated reports are
  outputs and are not committed.
- **No PDF:** a PDF library would be a heavy, fragile dependency. Print the
  HTML report to PDF from a browser instead.

### Contents

The report sections are:

1. Run Summary
2. Target, with provenance read from the PDB file's own header and the
   preparation log
3. Methodology, worded from what the run actually did
4. Screening Results: the top 3 and the full ranking
5. Docking Results
6. Molecular Descriptors
7. Unranked Compounds
8. Explainability Summary
9. Limitations
10. Responsible Use
11. Reproducibility

The JSON report has the top-level fields `metadata`, `target`,
`methodology`, `summary`, `compounds`, `ranking`, `explanations`,
`limitations`, `responsible_use` and `reproducibility`. NumPy values, paths
and datetimes are converted to plain JSON; NaN and infinity become `null`.

### Partial and failed runs

The report always describes what actually happened:

- **Vina unavailable:** descriptors are reported, and the docking and
  ranking sections say why they are empty.
- **No compounds ranked:** there is no Top Results section.
- **Some compounds invalid or failed:** they are listed with their status,
  the failed stage and the recorded error.
- **An empty run:** produces a valid empty report.

The pipeline reports the report's own status as `generated` (all formats
written), `partial` (some written) or `not_generated`. Writing a report can
never raise or discard the screening results already produced.

### Reproducibility metadata

Each report includes:

- **A traceable run ID:** UTC timestamp + target + a 10-character hash of
  every configured input that affects results (box, Vina settings, seed,
  ranking weights, docking scale, descriptor ranges, and the compound
  dataset and raw target checksums) + a 6-character hash of the run's
  outcomes. For example, `20261024-122417_6LU7_020ac65930-1f3a9c`. Two runs
  can share an ID only if their configuration *and* results are identical.
  No random UUIDs are used.
- **Software versions:** Python, RDKit, Meeko, gemmi, Streamlit and Vina.
- **The Git commit**, flagged if there are uncommitted changes.
- **The docking run folder.**
- **SHA-256 checksums** of the compound dataset, the raw target and the
  prepared receptor.

Prepared receptor files are written with LF line endings, the same way Git
stores them, so the receptor checksum matches a fresh checkout on any
platform.

### Limitations and responsible use

Every report carries the same fixed Limitations and Responsible Use
sections, and states that experimental validation is required. Reports
document computational descriptors, predicted docking results and
transparent prioritisation. They do not establish binding, efficacy,
toxicity, safety or clinical usefulness, and must not be used for clinical
decision-making. A test scans every format for unsupported claims.

## Project Architecture

| Folder | Role |
| ------ | ---- |
| `UI/` | Streamlit interface. `app.py` coordinates the app; `pages/` holds each page; `components/` holds reusable widgets; `utils/` holds display helpers. No scientific logic. |
| `backend/` | The screening pipeline. `screening.py` coordinates the stages: `compounds.py` (dataset loading), `descriptors.py`, `receptor.py` (target preparation and binding box), `ligands.py` (3D ligand preparation), `docking.py` (Vina), `ranking.py` (Screening Priority Score), `explainability.py` (rule-based explanations) and `report.py` (Markdown/JSON/HTML reports). No UI code. |
| `data/` | Input data: `compounds/` (compound libraries) and `targets/` (raw and prepared protein structures), each with provenance documentation. |
| `models/` | Trained ML models and their configuration. |
| `results/` | Generated docking runs, rankings and reports (ignored by Git). |
| `tests/` | Automated pytest tests: fast unit tests plus optional Vina integration tests. |
| `presentation/` | Hackathon slides, demo script and screenshots. |
| `config.py` | Single source of names, paths, limits, the binding box, docking settings, ranking weights and ranges, and logging. |

The UI only talks to the backend through functions such as
`run_screening()`, so each scientific stage can be developed and tested
independently.

## Folder Structure

```
Fairuz/
│
├── UI/
│   ├── app.py
│   ├── pages/
│   │   ├── home.py
│   │   ├── screening.py
│   │   └── results.py
│   ├── components/
│   │   ├── sidebar.py
│   │   ├── result_card.py
│   │   ├── explanation_panel.py
│   │   └── report_panel.py
│   └── utils/
│       └── helpers.py
│
├── backend/
│   ├── screening.py
│   ├── compounds.py
│   ├── descriptors.py
│   ├── receptor.py
│   ├── ligands.py
│   ├── docking.py
│   ├── ranking.py
│   ├── explainability.py
│   └── report.py
│
├── data/
│   ├── compounds/
│   │   └── compounds.csv
│   └── targets/
│       └── 6LU7/
│           ├── raw/6LU7.pdb
│           └── prepared/
│
├── models/
├── results/
├── tests/
│   ├── test_descriptors.py
│   ├── test_ranking.py
│   ├── test_compounds.py
│   ├── test_screening.py
│   ├── test_docking.py
│   ├── test_explainability.py
│   ├── test_report.py
│   ├── test_ui_helpers.py
│   └── conftest.py
│
├── presentation/
│   ├── DEMO_CHECKLIST.md
│   └── TECHNICAL_SUMMARY.md
├── config.py
├── pytest.ini
├── requirements.txt
├── README.md
├── .gitignore
└── LICENSE
```

Also present: `__init__.py` files (so folders import as packages), a
`README.md` in each data/output folder, and `.streamlit/config.toml` (theme).

## Technology Stack

**Current:** Python 3.11+, Streamlit, RDKit, Meeko, gemmi, AutoDock Vina
1.2.7, pytest

**Planned:** a supervised ML library (e.g. scikit-learn), only once a
suitable labelled dataset is added

## Installation

### 1. Python environment

```bash
git clone <repository-url>
cd Fairuz
python -m venv .venv
```

Activate it:

- Windows (PowerShell): `.venv\Scripts\Activate.ps1`
- macOS / Linux: `source .venv/bin/activate`

Install the Python dependencies (Streamlit, RDKit, Meeko, gemmi, SciPy,
pytest):

```bash
pip install -r requirements.txt
```

RDKit is installed from PyPI as `rdkit`. Do not use the old `rdkit-pypi`
package.

### 2. AutoDock Vina (separate download)

Vina is **not** installed by pip: its Python package has no prebuilt
Windows wheels and needs a C++/Boost toolchain to build. Fairouz instead
runs the official Vina executable.

1. Download the build for your system from the official release page:
   https://github.com/forlilab/AutoDock-Vina/releases/tag/v1.2.7
   (`vina_1.2.7_win.exe`, `vina_1.2.7_mac_aarch64`,
   `vina_1.2.7_mac_x86_64` or `vina_1.2.7_linux_x86_64`).
2. Make it findable, using **one** of these:
   - Save it into the virtual environment as `.venv/Scripts/vina.exe`
     (Windows) or `.venv/bin/vina` (macOS/Linux; then run `chmod +x`).
   - Put it on your PATH, named `vina`.
   - Set the environment variable `FAIROUZ_VINA` to its full path.
3. Check that it works: `vina --version` should print `AutoDock Vina v1.2.7`.

If Vina is missing, Fairouz still runs. Descriptors are calculated, and the
app explains that docking is unavailable and how to install Vina.

### 3. Prepared target (already included)

The prepared 6LU7 receptor is committed in `data/targets/6LU7/prepared/`.
To regenerate it from the raw RCSB file:

```bash
python -m backend.receptor
```

## Running Fairouz

From the project root:

```bash
streamlit run UI/app.py
```

The app opens at http://localhost:8501. Docking 10 compounds takes about
20–30 seconds.

## Testing

Fast unit tests (235 tests). These don't run Vina and take about 15 seconds:

```bash
python -m pytest
```

| Test file | Covers |
| --------- | ------ |
| `test_compounds.py` | CSV loading: missing columns, malformed rows, duplicate IDs, the real dataset |
| `test_descriptors.py` | RDKit descriptors against reference values, invalid SMILES, batches |
| `test_docking.py` | Ligand preparation, binding box, Vina command and output, setup errors, LF-stable receptor files |
| `test_ranking.py` | Screening Priority Score: scale, ranges, ties, unranked statuses, determinism |
| `test_explainability.py` | Explanations match the ranking exactly; wording; ties; failures |
| `test_report.py` | Report formats, schema, values, partial/failed runs, serialisation, run IDs, wording |
| `test_screening.py` | Pipeline stage statuses, timings, messages, degraded modes |
| `test_ui_helpers.py` | Display precision, no "−0.00", status icons |

Integration tests (5 tests). These run the real Vina executable,
regenerate the receptor from the raw PDB file, and run the full pipeline
twice to check that it is reproducible. They take about a minute and are
skipped if Vina is not installed:

```bash
python -m pytest -m integration
```

Everything:

```bash
python -m pytest -m ""
```

## Current Development Status

**Step 7 · Demo-ready prototype.** The full pipeline works end to end.

Implemented, in pipeline order:

1. **Compounds:** a PubChem-verified demo dataset, with checks for missing
   columns, malformed rows and duplicate IDs.
2. **Molecular analysis:** SMILES validation and five RDKit descriptors.
3. **Target:** SARS-CoV-2 main protease (PDB 6LU7), with documented
   provenance. Receptor preparation is reproducible and writes LF-stable
   files kept separate from the raw structure. The binding box is derived
   from the co-crystallised inhibitor N3.
4. **Ligands:** reproducible 3D preparation (fixed seeds, MMFF94, Meeko).
5. **Docking:** real AutoDock Vina, with failure handling per compound and
   each run saved with full metadata.
6. **Ranking:** a transparent Screening Priority Score (docking 70%,
   properties 30%), with explicit statuses for unranked compounds and
   deterministic tie-breaking.
7. **Explainability:** a rule-based breakdown of every score, taken exactly
   from the ranking, with factors, primary driver, comparisons and
   tie-break notes.
8. **Reports:** deterministic Markdown, JSON and HTML, with provenance,
   methods, limitations, a run ID, configuration and results hashes,
   versions and checksums.
9. **Pipeline status model:** every stage reports completed, partial,
   failed or skipped, with a one-line detail and its time. The outcome
   message is built from the final statuses.
10. **UI:**
    - **Home:** what Fairouz is, what is real and what it is not, the
      workflow, methodology and limitations.
    - **Screening:** readiness checks, live stage progress, stage status
      and timings.
    - **Results:** Overview, Explainability, Ranking Method, Details and
      Report tabs, with downloads.

Not yet implemented:

- Supervised ML ranking (deliberately not used; see "Why there is no
  machine-learning model yet")
- PDF reports (print the HTML report instead)
- Redocking validation with a non-covalent reference ligand (see
  `data/targets/6LU7/README.md`)

## Future Development

- Supervised ML ranking, only with a suitable labelled dataset (with
  model-appropriate explainability at that point)
- Redocking validation of the docking setup
- A larger, documented public compound library
- Comparing reports across screening runs (the JSON schema supports this)
- Final UX polish

## Responsible Use of Data

- **No patient data.** Fairouz works only with chemical and protein structure
  data. It does not require, collect or store patient data.
- **No personal information.** The app has no accounts and collects no
  personally identifiable information.
- **Public, documented data.** Compounds come from PubChem, and each row
  records its PubChem compound ID. The target structure comes from the RCSB
  Protein Data Bank (PDB 6LU7, CC0 public domain), with its citation,
  download date and checksum recorded in `data/targets/README.md`.
- **Inclusion is not endorsement.** A molecule's presence in the dataset does
  not imply any antiviral or therapeutic effect. The demo set exists only to
  test the software.
- **Using a viral target proves nothing by itself.** Docking a compound into
  the SARS-CoV-2 protease structure is a computational demonstration, not
  evidence of antiviral activity.
- **Problems are reported, not hidden.** Rows with missing values are skipped
  and listed, duplicate IDs stop loading, and molecules that fail parsing,
  3D generation or docking are marked failed with a reason. Nothing is
  silently accepted or filled in.
- **Descriptors are properties, not evidence.** Molecular weight, LogP and
  similar values describe a molecule. They do not show that it works as a
  drug.
- **Reproducibility.** Fixed random seeds, recorded tool versions, checksums
  of the input structures, and per-run metadata mean any docking result can
  be traced and repeated.
- **Small, unrepresentative demo set.** Ten hand-picked molecules do not
  represent the chemical space relevant to drug discovery.
- **Licensing.** Every external dataset must be checked for appropriate usage
  rights and licensing before it is added.
- **Bias.** Public chemical and biological datasets are often uneven: some
  chemical families and well-studied targets are over-represented. Rankings
  may favour compounds similar to those already well studied.
- **Training data limits.** Future ML models may inherit biases from their
  training datasets. A model can only be as reliable as the data it
  learns from.
- **Uncertainty.** Predictions should be reported with uncertainty where
  possible, so a ranking is never presented as more certain than it is.
- **Ranking is for research prioritisation only.** The Screening Priority
  Score is a transparent heuristic. Its weights and ranges are documented
  choices, not learned from data, and a high score is not evidence of
  efficacy, safety or clinical benefit.
- **No fake AI.** No model is trained on the 10 demo compounds, because the
  data cannot support one.
- **Not clinical evidence.** Computational results are predictions. They must
  not be treated as clinical evidence or as proof that a compound works.
- **Laboratory validation remains necessary** for every candidate.

## Interpreting Docking Scores

Fairouz reports each compound's **predicted docking affinity**: Vina's score
for its best pose, in kcal/mol.

- More negative values generally indicate a more favourable predicted pose
  *within the same docking setup* (same receptor, box and settings).
- The score is an approximation from a simplified scoring function. It is
  not a measured binding affinity.
- Comparing scores across studies, software, receptors or boxes can be
  misleading.
- Larger molecules tend to score more negatively simply because they make
  more contacts.
- One compound scoring better than another is not proof that it is
  biologically superior. Experimental testing remains essential.

## Molecular Docking Limitations

- Docking predicts likely *poses* and approximate interaction scores; it does
  not simulate real binding.
- The protein is treated as **rigid**, with simplified structural
  assumptions: no water molecules, a fixed histidine protonation state (HIE),
  and Gasteiger charges.
- A single crystal structure captures one conformation of the protein. Real
  proteins move, and 6LU7 is the shape the protease takes around one
  specific inhibitor.
- N3 is a covalent inhibitor; Vina models only non-covalent binding.
- Ligands are docked in the protonation and tautomer state written in their
  SMILES. Different states at physiological pH can change the results.
- The search box assumes compounds bind where N3 binds; other sites are not
  explored.
- Docking scores do not directly equal binding affinity.
- Docking scores do not prove biological activity.
- Docking scores do not establish safety.
- **Laboratory validation is necessary.**

## Limitations

- Fairouz is a hackathon proof-of-concept.
- The current dataset has only 10 molecules, chosen to test the software.
- Descriptor definitions vary between tools: RDKit's TPSA and H-bond counts
  can differ slightly from PubChem's for the same molecule.
- The screening pipeline is simplified compared with professional workflows.
- The docking setup has not yet been validated by redocking a reference
  ligand.
- The ranking inherits every docking limitation. It also adds its own
  assumptions (weights, the docking scale, property ranges) that have not
  been validated against experimental data.
- With 10 small demo molecules, most docking affinities are weaker than
  −6 kcal/mol, so the ranking spreads over a narrow band and
  property-component differences matter more than they would for drug-sized
  molecules.
- Any future ML model will depend on the quality and coverage of its
  training data.
- Nothing produced by Fairouz has been clinically validated.

## Disclaimer

Fairouz is a research and educational proof-of-concept. It is not intended for
clinical use, medical diagnosis, treatment decisions, or direct drug
development decisions without expert scientific review and experimental
validation.

## License

MIT. See [LICENSE](LICENSE).
