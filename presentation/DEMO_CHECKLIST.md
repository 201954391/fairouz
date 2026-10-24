# Fairouz demo checklist

Operational checklist for a live demo. Run through the "Before the demo"
section on the demo machine, ideally on the day.

## Before the demo

### 1. Environment

From the project root, with the virtual environment active:

```bash
python --version
```

Python 3.11 or newer.

```bash
python -c "import rdkit, meeko, gemmi, streamlit; print('imports ok')"
```

Prints `imports ok`.

### 2. AutoDock Vina

```bash
python -c "from backend.docking import find_vina_executable, vina_version; p = find_vina_executable(); print(p, vina_version(p))"
```

Expected: a path, then `AutoDock Vina v1.2.7`. If Vina is not found, see
README › Installation › 2. The app's sidebar also shows **AutoDock Vina found**
or **not found**.

### 3. Tests

```bash
python -m pytest
```

Expected: **235 passed** (about 15 s).

```bash
python -m pytest -m integration
```

Expected: **5 passed** (about 1 minute; runs real Vina).

### 4. Rehearsal run and backup report

Start the app, run one screening with the defaults, and confirm the expected
results below. This run also produces the backup report.

```bash
streamlit run UI/app.py
```

Expected results for the default run (target PDB 6LU7, 10 compounds):

| Check | Expected |
| ----- | -------- |
| Target | SARS-CoV-2 Main Protease (PDB 6LU7) |
| Compounds available | 10 |
| Pipeline stages | All 7 **Completed** |
| Runtime | About 20–30 s (docking is about 97% of it) |
| #1 | Ibuprofen: Screening Priority Score 54.2, predicted affinity −6.1 kcal/mol |
| #2 | Paracetamol: 43.9, −5.2 kcal/mol |
| #3 | Aspirin: 42.7, −5.1 kcal/mol |
| Tie | Benzene and Ethanol both 24.0; Benzene ranks higher (more negative affinity) |
| Report | Status **Generated**; Markdown, JSON and HTML downloads |

Then open the generated HTML report from `results/reports/` in a browser tab
and keep it open. That is the backup (see below). Note its run ID.

## Backup plan (if docking cannot run live)

Docking depends on the Vina executable and takes about 25 seconds, so have
a fallback that never misrepresents a result:

1. **Vina fails or is missing on the demo machine:** run the screening
   anyway. Fairouz still calculates descriptors, shows **Docking: Failed**
   and **Ranking: Skipped** with the reason, and generates a report saying
   so. This demonstrates the failure handling honestly.
2. **Then show the backup report** you opened in the browser, and say so
   explicitly: *"This is a report Fairouz generated earlier on this machine,
   run ID …, from a real Vina run."* Its run ID, timestamp, Git commit and
   checksums identify it.
3. **Never** present the backup as a live run, and never edit a report by
   hand.

## Demo path (about 3 minutes)

| Time | Show | Say (in your own words) |
| ---- | ---- | ----------------------- |
| 0:00–0:20 | **Home** | Testing compounds in the lab is slow and expensive. Fairouz prioritises which ones are worth investigating, transparently. |
| 0:20–0:45 | Home: target line, workflow, "What is real" | A real X-ray structure (PDB 6LU7), 10 public molecules verified against PubChem, real RDKit and AutoDock Vina. |
| 0:45–1:15 | **Screening** → Start screening | Point at the live stage messages: descriptors, 3D ligands, Vina compound by compound, scores, explanations, report. Then the 7 stage statuses and the runtime. |
| 1:15–1:40 | **Results → Overview** | Top 3 cards: Screening Priority Score and predicted docking affinity. One disclaimer: computational screening output, not evidence of efficacy. |
| 1:40–2:10 | **Explainability** (pick Ibuprofen, then Benzene) | 24.2 docking points + 30.0 property points = 54.2, straight from the formula. The tie-break note for Benzene and Ethanol. No language model. |
| 2:10–2:30 | **Ranking Method**, then **Report** | The formula and ranges, no hidden weights. Download the HTML report: run ID, provenance, limitations, checksums. |
| 2:30–3:00 | Close | No ML on purpose: 10 molecules give no labelled training data, so Fairouz uses transparent deterministic prioritisation. 235 unit tests plus real-Vina integration tests. Every result needs experimental validation. |

## Key tabs to have ready

1. Home
2. Screening
3. Results → Overview
4. Results → Explainability
5. Results → Ranking Method
6. Results → Report

## Closing message

> Fairouz shows that a virtual-screening workflow can be both real and
> transparent: real structures and real docking, with every score explained
> and every run reproducible. It prioritises compounds for further
> investigation. It does not show that any compound works, and every result
> needs experimental validation.
