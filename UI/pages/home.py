"""Home page: what Fairouz is, what is real, and how to start."""

import streamlit as st

import config
from UI.utils.helpers import (
    library_size,
    render_disclaimer,
    render_limitations,
    render_workflow,
    target_provenance,
)

# Set by UI/app.py so the call-to-action can link to the Screening page.
screening_page: st.Page | None = None

NO_ML_STATEMENT = (
    "Fairouz currently uses transparent deterministic prioritisation because the "
    "available compound set does not provide an appropriate labelled training dataset."
)


def _target_line() -> str:
    facts = target_provenance()
    method = (facts.get("method") or "").capitalize()
    resolution = facts.get("resolution_angstrom")
    detail = ", ".join(part for part in (method, f"{resolution:g} Å" if resolution else "") if part)
    return f"**Target:** {config.TARGET_LABEL}" + (f" · {detail}" if detail else "")


def _render_technical_details() -> None:
    """Compact implementation details for judges; collapsed by default."""
    weak, strong = config.DOCKING_AFFINITY_WEAK, config.DOCKING_AFFINITY_STRONG
    docking_pct, property_pct = round(100 * config.RANKING_DOCKING_WEIGHT), round(100 * config.RANKING_PROPERTY_WEIGHT)
    with st.expander("About & methodology (technical details)"):
        st.markdown(
            f"""
- **Target:** PDB {config.TARGET_ID} from the {config.TARGET_SOURCE}; protease chain {config.RECEPTOR_CHAIN} kept,
  the co-crystallised inhibitor N3 and waters removed, polar hydrogens and charges added with Meeko.
- **Binding box:** derived from N3's coordinates with {config.BINDING_BOX_PADDING:g} Å padding; centre
  {config.BINDING_BOX_CENTER} Å, size {config.BINDING_BOX_SIZE} Å.
- **Compounds:** public molecules whose SMILES strings were verified against PubChem; each row records its PubChem ID.
- **Molecular analysis:** RDKit calculates molecular weight, Wildman–Crippen LogP, TPSA and hydrogen-bond donors/acceptors.
- **Ligands:** RDKit ETKDGv3 3D conformer (seed {config.RANDOM_SEED}), MMFF94 optimisation, Meeko PDBQT conversion.
- **Docking:** the AutoDock Vina executable; exhaustiveness {config.VINA_EXHAUSTIVENESS},
  {config.VINA_NUM_POSES} poses, seed {config.RANDOM_SEED}.
- **Ranking:** Screening Priority Score = docking points (up to {docking_pct}) + property points (up to {property_pct}).
  Docking uses a fixed scale ({weak:g} kcal/mol → 0, {strong:g} kcal/mol → maximum); properties use documented
  Lipinski/Ghose/Veber-style ranges.
- **Explainability:** fixed rules applied to the ranking's own numbers; no language model.
- **Reports:** Markdown, JSON and HTML, each with a traceable run ID, configuration hash, software versions and checksums.
- **Provenance:** `data/compounds/README.md` and `data/targets/README.md`.

{NO_ML_STATEMENT}
"""
        )


def render() -> None:
    """Render the home page."""
    st.markdown(
        f'<p class="fairouz-wordmark">{config.PROJECT_NAME.upper()}</p>'
        f'<p class="fairouz-subtitle">{config.PROJECT_SUBTITLE}</p>',
        unsafe_allow_html=True,
    )
    st.write(config.PROJECT_DESCRIPTION)
    size = library_size()
    st.markdown(_target_line() + (f" · **Compounds:** {size} public molecules" if size else ""))
    if screening_page is not None and st.button("Start a screening", type="primary", icon=":material/play_arrow:"):
        st.switch_page(screening_page)

    does, real, not_claimed = st.columns(3)
    with does, st.container(border=True):
        st.markdown("**What it does**")
        st.write(
            "Docks each compound into the viral protein, scores it transparently, and ranks the set "
            "so researchers can choose which compounds to investigate further."
        )
    with real, st.container(border=True):
        st.markdown("**What is real**")
        st.write(
            "A real X-ray protein structure, RDKit descriptor calculations and AutoDock Vina docking. "
            "Every number shown is computed by the pipeline; nothing is simulated."
        )
    with not_claimed, st.container(border=True):
        st.markdown("**What it is not**")
        st.write(
            "Not evidence that any compound works. Docking scores are computational estimates, and "
            "every result needs experimental validation."
        )

    st.subheader("Workflow")
    render_workflow(config.PIPELINE_STAGES)

    _render_technical_details()
    with st.expander("Limitations"):
        render_limitations()

    render_disclaimer()
