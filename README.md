# Fairouz

**AI-Assisted Virtual Drug Screening** · research proof-of-concept

> Fairouz is a research and educational proof-of-concept. It is not intended
> for clinical use, medical diagnosis, treatment decisions, or direct drug
> development decisions without expert scientific review and experimental
> validation.

## Overview

Fairouz is a virtual drug screening platform. A researcher selects a viral
protein target, and Fairouz will screen a library of chemical compounds
against it, rank them with AI and explain why each one ranked highly.

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
learning, it produces a ranked shortlist of compounds together with the
reasons behind each ranking and the uncertainty of each prediction.
Researchers can then focus experimental effort on a smaller, better-justified
set of candidates. The aim is to reduce unnecessary screening effort, not to
replace experiments.

## How It Will Work

```
Viral Protein Target
  → Compound Library
  → Molecular Descriptors
  → Molecular Docking
  → AI Ranking
  → Top Candidates
  → Research Report
  → Laboratory Validation
```

| Stage | Planned approach |
| ----- | ---------------- |
| Target | Prepared public protein structure (e.g. from the Protein Data Bank) |
| Compound library | Curated public compound dataset with documented provenance |
| Molecular descriptors | RDKit: molecular weight, LogP, TPSA, H-bond donors/acceptors |
| Molecular docking | AutoDock Vina: predicted binding affinity (kcal/mol) |
| AI ranking | Simple, explainable ML model combining descriptors and docking scores |
| Report | Research-style summary with limitations and uncertainty |
| Laboratory validation | Carried out by researchers; outside the software |

## Project Architecture

| Folder | Role |
| ------ | ---- |
| `UI/` | Streamlit interface. `app.py` coordinates the app; `pages/` holds each page; `components/` holds reusable widgets; `utils/` holds display helpers. No scientific logic. |
| `backend/` | The screening pipeline. `screening.py` coordinates the stages in `descriptors.py`, `docking.py`, `ranking.py` and `report.py`. No UI code. |
| `data/` | Input data: `compounds/` (compound libraries) and `targets/` (prepared protein structures), each with a provenance table. |
| `models/` | Trained ML models and their configuration. |
| `results/` | Generated screening results, rankings and reports (ignored by Git). |
| `tests/` | Automated pytest tests for the backend. |
| `presentation/` | Hackathon slides, demo script and screenshots. |
| `config.py` | Single source of names, paths, limits and logging settings. |

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
│   ├── docking.py
│   ├── descriptors.py
│   ├── ranking.py
│   └── report.py
│
├── data/
│   ├── compounds/
│   └── targets/
│
├── models/
├── results/
├── tests/
│   ├── test_descriptors.py
│   └── test_ranking.py
│
├── presentation/
├── config.py
├── requirements.txt
├── README.md
├── .gitignore
└── LICENSE
```

Also present: `__init__.py` files (so folders import as packages), a
`README.md` in each data/output folder, and `.streamlit/config.toml` (theme).

## Technology Stack

**Current:** Python 3.11+, Streamlit, pytest

**Planned:** RDKit, AutoDock Vina, a simple ML library (e.g. scikit-learn)

## Installation

```bash
git clone <repository-url>
cd Fairuz
python -m venv .venv
```

Activate the virtual environment:

- Windows (PowerShell): `.venv\Scripts\Activate.ps1`
- macOS / Linux: `source .venv/bin/activate`

Install dependencies:

```bash
pip install -r requirements.txt
```

## Running Fairouz

From the project root:

```bash
streamlit run UI/app.py
```

The app opens at http://localhost:8501.

## Testing

From the project root:

```bash
python -m pytest
```

## Current Development Status

**Step 1 · Foundation & Interface.** This step builds the project structure,
configuration, logging and a working Streamlit interface. The backend modules
are placeholders that report "not implemented". No scientific calculations
are performed and no results, scores or compounds are generated.

## Future Development

- Real, documented public compound data
- RDKit molecular descriptor calculation
- Viral protein target preparation
- Molecular docking with AutoDock Vina
- ML-based ranking
- Explainability (why each compound ranked highly)
- Research report generation
- Final UX polish

## Responsible Use of Data

- **No patient data.** Fairouz works only with chemical and protein structure
  data. It does not require, collect or store patient or personal data.
- **Public data.** The project intends to use publicly available scientific
  and chemical datasets.
- **Licensing.** Every external dataset must be checked for appropriate usage
  rights and licensing before it is added.
- **Provenance.** The source, version, licence and retrieval date of every
  dataset and target structure is recorded in the provenance tables in
  `data/compounds/README.md` and `data/targets/README.md`.
- **Quality and completeness.** Dataset quality, missing values and
  inconsistent measurements must be checked and documented.
- **Bias.** Public chemical and biological datasets are often uneven: some
  chemical families and well-studied targets are over-represented. Rankings
  may favour compounds similar to those already well studied.
- **Training data limits.** A model can only be as reliable as the data it
  learns from. Small or unrepresentative training data will limit how well
  predictions generalise to new compounds.
- **Uncertainty.** Predictions should be reported with uncertainty where
  possible, so a ranking is never presented as more certain than it is.
- **Not clinical evidence.** Computational results are predictions. They must
  not be treated as clinical evidence or as proof that a compound works.
- **Laboratory validation remains necessary** for every candidate.

## Limitations

- Fairouz is a hackathon proof-of-concept.
- It will start with a small dataset.
- The screening pipeline is simplified compared with professional workflows.
- Docking scores are estimates of binding, not measurements.
- ML performance depends on the quality and coverage of the training data.
- Nothing produced by Fairouz has been clinically validated.

## Disclaimer

Fairouz is a research and educational proof-of-concept. It is not intended for
clinical use, medical diagnosis, treatment decisions, or direct drug
development decisions without expert scientific review and experimental
validation.

## License

MIT. See [LICENSE](LICENSE).
