# Fairouz

**AI-Assisted Virtual Drug Screening** · research proof-of-concept

> Fairouz is a research and educational proof-of-concept. It is not intended
> for clinical use, medical diagnosis, treatment decisions, or direct drug
> development decisions without expert scientific review and experimental
> validation.

## Overview

Fairouz is a virtual drug screening platform. A researcher selects a viral
protein target, and Fairouz screens a library of chemical compounds against
it. The final version will also rank the compounds with AI and explain why
each one ranked highly.

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
By combining molecular property analysis, molecular docking and machine
learning, it will produce a ranked shortlist of compounds together with the
reasons behind each ranking and the uncertainty of each prediction.
Researchers can then focus experimental effort on a smaller, better-justified
set of candidates. The aim is to reduce unnecessary screening effort, not to
replace experiments.

## How It Works

```
Viral Protein Target
  → Compound Library
  → Molecular Descriptors
  → 3D Ligand Preparation
  → Molecular Docking
  → AI Ranking
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
| AI ranking | Simple, explainable ML model combining descriptors and docking scores | Planned |
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
5. **Storage**: each run gets a new timestamped folder in `results/docking/`
   containing `docking_results.csv`, `run_metadata.json` (all settings,
   checksums and tool versions), and the ligand and pose files.

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

## Project Architecture

| Folder | Role |
| ------ | ---- |
| `UI/` | Streamlit interface. `app.py` coordinates the app; `pages/` holds each page; `components/` holds reusable widgets; `utils/` holds display helpers. No scientific logic. |
| `backend/` | The screening pipeline. `screening.py` coordinates the stages: `compounds.py` (dataset loading), `descriptors.py`, `receptor.py` (target preparation and binding box), `ligands.py` (3D ligand preparation), `docking.py` (Vina), `ranking.py` and `report.py`. No UI code. |
| `data/` | Input data: `compounds/` (compound libraries) and `targets/` (raw and prepared protein structures), each with provenance documentation. |
| `models/` | Trained ML models and their configuration. |
| `results/` | Generated docking runs, rankings and reports (ignored by Git). |
| `tests/` | Automated pytest tests: fast unit tests plus optional Vina integration tests. |
| `presentation/` | Hackathon slides, demo script and screenshots. |
| `config.py` | Single source of names, paths, limits, the binding box, docking settings and logging. |

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
│   │   ├── screening.py
│   │   └── results.py
│   ├── components/
│   │   ├── sidebar.py
│   │   └── result_card.py
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
│   └── test_docking.py
│
├── presentation/
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

**Planned:** a simple ML library (e.g. scikit-learn)

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

Fast unit tests. These don't run Vina and take about 2 seconds:

```bash
python -m pytest
```

Docking integration tests. These run the real Vina executable, take about
15 seconds, and are skipped if Vina is not installed:

```bash
python -m pytest -m integration
```

Everything:

```bash
python -m pytest -m ""
```

## Current Development Status

**Step 3 · Molecular Docking.**

Implemented:

- Compound dataset loading and SMILES validation (Step 2)
- RDKit molecular descriptors and their display (Step 2)
- Real viral target: SARS-CoV-2 main protease, PDB 6LU7, with documented
  provenance
- Reproducible receptor preparation, with raw and prepared files kept separate
- Binding box derived from the co-crystallised inhibitor N3
- Reproducible 3D ligand preparation (fixed seeds)
- Real AutoDock Vina docking with per-compound failure handling
- Docking results saved per run with full reproducibility metadata
- Docking results table, compound cards and score explanations in the UI

Not yet implemented:

- AI ranking
- Report generation
- Redocking validation with a non-covalent reference ligand (see
  `data/targets/6LU7/README.md`)

## Future Development

- AI/ML candidate ranking using descriptor and docking data
- Explainability (why each compound ranked highly)
- Redocking validation of the docking setup
- A larger, documented public compound library
- Research report generation
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
- ML performance will depend on the quality and coverage of the training
  data.
- Nothing produced by Fairouz has been clinically validated.

## Disclaimer

Fairouz is a research and educational proof-of-concept. It is not intended for
clinical use, medical diagnosis, treatment decisions, or direct drug
development decisions without expert scientific review and experimental
validation.

## License

MIT. See [LICENSE](LICENSE).
