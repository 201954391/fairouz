"""Results page: docking results, molecular properties and future sections."""

from typing import Any

import streamlit as st

import config
from UI.components.result_card import render_result_card
from UI.pages.screening import LAST_SCREENING_KEY
from UI.utils.helpers import (
    AFFINITY_EXPLANATION,
    AFFINITY_HELP,
    AFFINITY_LABEL,
    descriptor_label,
    format_status,
    render_disclaimer,
)

AFFINITY_COLUMN = f"{AFFINITY_LABEL} (kcal/mol)"
SORT_BY_AFFINITY = "Predicted affinity (most negative first)"
SORT_BY_DATASET = "Dataset order"

SCORE_GUIDANCE = """
- **More negative is more favourable** — but only relative to other compounds
  docked with *this* receptor, box and settings.
- **It is an approximation.** Vina's scoring function is a simplified model of
  binding; the number is not a measured binding affinity.
- **Do not compare across studies.** Scores from other software, receptors or
  boxes are not on the same footing.
- **Size matters.** Larger molecules tend to score more negatively simply
  because they make more contacts.
- **It is not evidence of activity or safety.** Experimental testing remains
  essential.
"""

# Integer-valued descriptors are shown without decimals.
COUNT_DESCRIPTORS = {"h_bond_donors", "h_bond_acceptors"}

# Sections whose backend stage is not built yet: tab title -> (stage key, description).
PENDING_SECTIONS = {
    "Top Ranked Compounds": (
        "ai_ranking",
        "The compounds Fairouz predicts are most worth investigating further, "
        "with the main reasons each one ranked highly.",
    ),
    "AI Ranking": (
        "ai_ranking",
        "How the ranking model combined the features, including the "
        "uncertainty of each prediction.",
    ),
    "Final Report": (
        "report",
        "A downloadable research summary covering the target, data sources, "
        "top candidates, limitations and next steps for laboratory validation.",
    ),
}
IMPLEMENTED_SECTIONS = ("Docking Results", "Molecular Properties")
TAB_ORDER = ("Top Ranked Compounds", "Docking Results", "Molecular Properties", "AI Ranking", "Final Report")


def _descriptor_rows(compounds: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Build one table row per compound, using only RDKit-calculated values."""
    rows = []
    for compound in compounds:
        row = {"Name": compound["name"]}
        for key in config.DESCRIPTORS:
            row[descriptor_label(key)] = compound[key]
        row["Status"] = "Calculated" if compound["valid"] else f"Failed: {compound['error']}"
        rows.append(row)
    return rows


def _descriptor_column_config() -> dict[str, Any]:
    """Column formats and tooltips for the descriptor table."""
    columns: dict[str, Any] = {
        "Status": st.column_config.TextColumn(
            "Status", help="Whether RDKit could calculate the descriptors."
        ),
    }
    for key, descriptor in config.DESCRIPTORS.items():
        columns[descriptor_label(key)] = st.column_config.NumberColumn(
            descriptor_label(key),
            help=descriptor["help"],
            format="%d" if key in COUNT_DESCRIPTORS else "%.2f",
        )
    return columns


def _docking_status(compound: dict[str, Any]) -> str:
    """Short, honest docking status for one compound."""
    docking = compound.get("docking")
    if docking is None:
        return "Not docked: invalid molecule" if not compound["valid"] else "Not docked"
    if docking["docking_success"]:
        return "Docked"
    return f"Failed: {docking['error']}"


def _docking_rows(compounds: list[dict[str, Any]], sort_by_affinity: bool) -> list[dict[str, Any]]:
    """Build the docking table using only affinities returned by Vina."""
    rows = [
        {
            "Compound": compound["name"],
            AFFINITY_COLUMN: (compound.get("docking") or {}).get("best_affinity_kcal_mol"),
            "Docking Status": _docking_status(compound),
            "Poses": (compound.get("docking") or {}).get("poses_generated"),
        }
        for compound in compounds
    ]
    if sort_by_affinity:
        # Docked compounds first, most negative affinity first; the rest keep order.
        rows.sort(key=lambda row: (row[AFFINITY_COLUMN] is None, row[AFFINITY_COLUMN] or 0.0))
    return rows


def _render_docking_setup(docking_info: dict[str, Any]) -> None:
    """Show the settings needed to reproduce the docking run."""
    metadata = docking_info["metadata"]
    box, vina = metadata["binding_box"], metadata["vina"]
    with st.expander("Docking setup (for reproducibility)"):
        st.markdown(
            f"- **Target:** {config.TARGET_LABEL}, receptor chain {config.RECEPTOR_CHAIN}\n"
            f"- **Engine:** {vina['version']}, exhaustiveness {vina['exhaustiveness']}, "
            f"{vina['num_modes']} poses, seed {vina['seed']}\n"
            f"- **Binding box (Å):** centre {tuple(box['center'])}, size {tuple(box['size'])}, "
            "derived from the co-crystallised inhibitor N3\n"
            f"- **Ligands:** {metadata['ligand_preparation']['conformer']}, "
            f"{metadata['ligand_preparation']['force_field']}\n"
            f"- **Saved to:** `{docking_info['run_dir']}`"
        )


def _render_docking_results(last_screening: dict[str, Any]) -> None:
    """Show the computational docking results table."""
    docking_info = last_screening["docking"]
    if docking_info["error"]:
        st.warning(f"Docking could not run: {docking_info['error']}", icon=":material/warning:")
        return

    st.markdown("**Computational Docking Results**")
    st.caption(AFFINITY_EXPLANATION)
    sort_choice = st.radio(
        "Sort by", (SORT_BY_AFFINITY, SORT_BY_DATASET), horizontal=True, label_visibility="collapsed"
    )
    st.dataframe(
        _docking_rows(last_screening["compounds"], sort_choice == SORT_BY_AFFINITY),
        column_config={
            AFFINITY_COLUMN: st.column_config.NumberColumn(AFFINITY_COLUMN, help=AFFINITY_HELP, format="%.1f"),
            "Docking Status": st.column_config.TextColumn("Docking Status"),
            "Poses": st.column_config.NumberColumn("Poses", help="Poses returned by Vina.", format="%d"),
        },
        hide_index=True,
        width="stretch",
    )
    with st.expander("How to read this score"):
        st.markdown(SCORE_GUIDANCE)
    _render_docking_setup(docking_info)


def _render_molecular_properties(compounds: list[dict[str, Any]]) -> None:
    """Show the descriptor table and one card per compound."""
    st.write(
        "Physical and chemical properties calculated by RDKit. Hover over a "
        "column heading for a short explanation."
    )
    st.caption(
        "These are properties of each molecule, not evidence that it is "
        "effective against any virus."
    )
    st.dataframe(
        _descriptor_rows(compounds),
        column_config=_descriptor_column_config(),
        hide_index=True,
        width="stretch",
    )

    with st.expander(f"Compound cards ({len(compounds)})"):
        columns = st.columns(2)
        for index, compound in enumerate(compounds):
            with columns[index % 2]:
                render_result_card(
                    compound_name=compound["name"],
                    properties={key: compound[key] for key in config.DESCRIPTORS},
                    smiles=compound["smiles"],
                    error=compound["error"],
                    docking=compound.get("docking"),
                )


def _render_pending_section(stage_status: str, description: str) -> None:
    """Show a section whose backend stage is not built yet."""
    st.write(description)
    st.info(f"Status: **{format_status(stage_status)}**. No values are shown until this stage is built.")


def render() -> None:
    """Render the results page."""
    st.title("Results")
    last_screening = st.session_state.get(LAST_SCREENING_KEY)

    if last_screening is None:
        st.info(
            "No screening has been run yet. Start a screening to calculate "
            "molecular properties.",
            icon=":material/hourglass_empty:",
        )
        stages: dict[str, str] = {}
    else:
        request, summary = last_screening["request"], last_screening["summary"]
        st.write(
            f"Latest run: **{request['target']}** · {summary['processed']} compounds "
            f"processed · {summary['valid']} analysed · {summary['docked']} docked · "
            f"{summary['invalid'] + summary['docking_failed']} failed."
        )
        stages = last_screening["stages"]

    for tab, title in zip(st.tabs(TAB_ORDER), TAB_ORDER):
        with tab:
            st.markdown(f"#### {title}")
            if title in IMPLEMENTED_SECTIONS:
                if last_screening is None:
                    st.caption("No data available yet.")
                elif title == "Docking Results":
                    _render_docking_results(last_screening)
                else:
                    _render_molecular_properties(last_screening["compounds"])
            else:
                stage, description = PENDING_SECTIONS[title]
                _render_pending_section(stages.get(stage, config.STATUS_NOT_IMPLEMENTED), description)

    render_disclaimer()
