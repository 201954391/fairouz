"""Results page: overview, explainability, ranking method, details and report."""

from typing import Any

import streamlit as st

import config
from UI.components.explanation_panel import overview_chart, render_explanation
from UI.components.report_panel import render_report
from UI.components.result_card import render_result_card, render_top_card
from UI.pages.home import NO_ML_STATEMENT
from UI.pages.screening import LAST_SCREENING_KEY
from UI.utils.helpers import (
    AFFINITY_EXPLANATION,
    AFFINITY_HELP,
    AFFINITY_LABEL,
    PRIORITY_EXPLANATION,
    PRIORITY_HELP,
    PRIORITY_LABEL,
    descriptor_label,
    descriptor_number_format,
    format_affinity,
    format_priority,
    format_seconds,
    render_limitations,
    render_pipeline_status,
    round_for_display,
)

TOP_CARD_COUNT = 3
TAB_ORDER = ("Overview", "Explainability", "Ranking Method", "Details", "Report")

AFFINITY_COLUMN = f"{AFFINITY_LABEL} (kcal/mol)"
AFFINITY_FORMAT = f"%.{config.AFFINITY_DECIMALS}f"
SCORE_FORMAT = f"%.{config.SCORE_DECIMALS}f"
SORT_BY_AFFINITY = "Predicted affinity (most negative first)"
SORT_BY_DATASET = "Dataset order"

SCORE_GUIDANCE = """
- **More negative is more favourable**, but only relative to other compounds
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


# --- Shared row builders -----------------------------------------------------------

def _descriptor_values(compound: dict[str, Any]) -> dict[str, Any]:
    """Descriptor table cells, rounded for display (no "-0.00")."""
    return {
        descriptor_label(key): round_for_display(compound[key], descriptor["decimals"])
        for key, descriptor in config.DESCRIPTORS.items()
    }


def _descriptor_column_config() -> dict[str, Any]:
    return {
        descriptor_label(key): st.column_config.NumberColumn(
            descriptor_label(key), help=descriptor["help"], format=descriptor_number_format(key)
        )
        for key, descriptor in config.DESCRIPTORS.items()
    }


def _ranked(last_screening: dict[str, Any]) -> list[dict[str, Any]]:
    return [entry for entry in last_screening["ranking"]["ranked_compounds"] if entry["rank"] is not None]


def _compounds_by_id(last_screening: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {compound["compound_id"]: compound for compound in last_screening["compounds"]}


# --- Overview ----------------------------------------------------------------------

def _render_headline(last_screening: dict[str, Any], ranked: list[dict[str, Any]]) -> None:
    summary = last_screening["summary"]
    st.markdown(f"**Target:** {last_screening['request']['target']}")
    processed, docked, ranked_count = st.columns(3)
    processed.metric("Processed", summary["processed"])
    docked.metric("Docked", summary["docked"], help="Compounds AutoDock Vina docked successfully.")
    ranked_count.metric("Ranked", summary["ranked"], help="Compounds given a Screening Priority Score.")
    if ranked:
        best = ranked[0]
        st.markdown(
            f"**Top-ranked compound: {best['name']}** · {PRIORITY_LABEL} "
            f"**{format_priority(best['screening_priority_score'])}** · predicted docking affinity "
            f"**{format_affinity(best['docking_affinity'])}**",
            help=f"{PRIORITY_HELP}\n\n{AFFINITY_HELP}",
        )


def _ranked_rows(ranked: list[dict[str, Any]], compounds_by_id: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "Rank": entry["rank"],
            "Compound": entry["name"],
            PRIORITY_LABEL: round_for_display(entry["screening_priority_score"], config.SCORE_DECIMALS),
            AFFINITY_COLUMN: round_for_display(entry["docking_affinity"], config.AFFINITY_DECIMALS),
            **_descriptor_values(compounds_by_id[entry["compound_id"]]),
            "Why": entry["explanation"],
        }
        for entry in ranked
    ]


def _ranked_column_config() -> dict[str, Any]:
    return {
        "Rank": st.column_config.NumberColumn("Rank", format="%d"),
        PRIORITY_LABEL: st.column_config.ProgressColumn(
            PRIORITY_LABEL, help=PRIORITY_HELP, min_value=0, max_value=100, format=SCORE_FORMAT
        ),
        AFFINITY_COLUMN: st.column_config.NumberColumn(AFFINITY_COLUMN, help=AFFINITY_HELP, format=AFFINITY_FORMAT),
        **_descriptor_column_config(),
        "Why": st.column_config.TextColumn("Why", help="Rule-based summary of what drove the score."),
    }


def _render_overview(last_screening: dict[str, Any]) -> None:
    ranked = _ranked(last_screening)
    compounds_by_id = _compounds_by_id(last_screening)
    unranked = [entry for entry in last_screening["ranking"]["ranked_compounds"] if entry["rank"] is None]

    _render_headline(last_screening, ranked)
    st.info(PRIORITY_EXPLANATION, icon=":material/info:")
    render_pipeline_status(last_screening["pipeline"])

    if not ranked:
        st.warning(
            "No compound was ranked: ranking needs successful docking results. "
            "Molecular descriptors are still available under Details.",
            icon=":material/warning:",
        )
    else:
        st.markdown(f"**Highest screening priority** (top {min(TOP_CARD_COUNT, len(ranked))})")
        for column, entry in zip(st.columns(TOP_CARD_COUNT), ranked[:TOP_CARD_COUNT]):
            with column:
                render_top_card(entry)
        st.markdown("**Full ranking**")
        st.dataframe(
            _ranked_rows(ranked, compounds_by_id), column_config=_ranked_column_config(), hide_index=True, width="stretch"
        )

    if unranked:
        st.markdown("**Not ranked**")
        st.dataframe(
            [{"Compound": entry["name"], "Status": entry["ranking_status"], "Reason": entry["explanation"]} for entry in unranked],
            hide_index=True,
            width="stretch",
        )

    st.caption(
        f"{last_screening['summary']['processed']} compounds processed in "
        f"{format_seconds(last_screening['timings']['total'])}. Stage timings are on the Screening page."
    )
    with st.expander("Limitations"):
        render_limitations()


# --- Explainability -----------------------------------------------------------------

def _explanation_label(explanation: dict[str, Any]) -> str:
    if explanation["ranked"]:
        return f"#{explanation['rank']} · {explanation['name']} ({explanation['score_display']})"
    return f"Not ranked · {explanation['name']}"


def _render_explainability(last_screening: dict[str, Any]) -> None:
    explanations = last_screening["explanations"]
    st.write(
        "Every number below comes directly from the Screening Priority Score calculation. "
        "Explanations are generated by fixed rules, not by a language model, so the same "
        "ranking always gives the same explanation."
    )
    if not explanations:
        st.caption("No compounds to explain.")
        return
    if any(explanation["ranked"] for explanation in explanations):
        st.markdown("**Ranking at a glance**")
        st.altair_chart(overview_chart(explanations), width="stretch")

    labels = [_explanation_label(explanation) for explanation in explanations]
    # A run-specific key, so a new run never shows a stale selection.
    run_id = (last_screening.get("report") or {}).get("run_id") or "latest"
    choice = st.selectbox("Explain compound", labels, key=f"explain_compound_{run_id}")
    render_explanation(explanations[labels.index(choice)])


# --- Ranking method ------------------------------------------------------------------

def _render_ranking_method(last_screening: dict[str, Any]) -> None:
    method = last_screening["ranking"]["method"]
    scale = method["docking_scale_kcal_mol"]
    docking_max, property_max = round(100 * method["docking_weight"]), round(100 * method["property_weight"])
    st.markdown(
        f"**{PRIORITY_LABEL}** = docking points + property points (0 to 100).\n\n"
        f"- **Docking points (up to {docking_max})**: the predicted Vina affinity on a fixed scale, where "
        f"{scale['zero_at']:g} kcal/mol or weaker earns 0 and {scale['full_at']:g} kcal/mol or stronger earns "
        f"all {docking_max}; linear in between. Docking is the strongest signal.\n"
        f"- **Property points (up to {property_max})**: the average of five descriptor sub-scores × {property_max}. "
        "Each sub-score is 1 inside its configured range and falls linearly to 0 outside it.\n"
        "- **Ties** are ordered by more negative docking affinity, then compound ID."
    )
    st.dataframe(
        [
            {
                "Property": config.DESCRIPTORS[key]["label"],
                # A range whose zero and full limits coincide has no lower penalty.
                "Scores 0 at or below": None if zero_below == full_from else zero_below,
                "Full score from": full_from,
                "Full score to": full_to,
                "Scores 0 at or above": zero_above,
            }
            for key, (zero_below, full_from, full_to, zero_above) in method["property_ranges"].items()
        ],
        hide_index=True,
        width="stretch",
    )
    st.caption(
        "Full-score ranges follow common rules of thumb (Lipinski, Ghose, Veber) for "
        "typical oral drugs. They do not show that a molecule is suitable or effective."
    )

    ranked = _ranked(last_screening)
    if ranked:
        st.markdown("**Score breakdown per compound**")
        st.dataframe(
            [
                {
                    "Rank": entry["rank"],
                    "Compound": entry["name"],
                    "Docking points": round_for_display(entry["docking_points"], 1),
                    "Property points": round_for_display(entry["property_points"], 1),
                    PRIORITY_LABEL: round_for_display(entry["screening_priority_score"], config.SCORE_DECIMALS),
                    **{
                        f"{config.DESCRIPTORS[key]['label']} sub-score": round_for_display(score, 2)
                        for key, score in entry["property_subscores"].items()
                    },
                }
                for entry in ranked
            ],
            column_config={"Rank": st.column_config.NumberColumn("Rank", format="%d")},
            hide_index=True,
            width="stretch",
        )

    with st.expander("Why is there no machine-learning model?"):
        st.markdown(
            f"{NO_ML_STATEMENT} A supervised model trained on a handful of molecules would memorise noise, "
            "so its output would be misleading. One may be added later, but only with a documented, licensed, "
            "labelled dataset (for example published Mpro activity measurements) and a held-out evaluation."
        )


# --- Details: docking and molecular properties ---------------------------------------

def _docking_status(compound: dict[str, Any]) -> str:
    docking = compound.get("docking")
    if docking is None:
        return "Not docked: invalid molecule" if not compound["valid"] else "Not docked"
    if docking["docking_success"]:
        return "Docked"
    return f"Failed: {docking['error']}"


def _docking_rows(compounds: list[dict[str, Any]], sort_by_affinity: bool) -> list[dict[str, Any]]:
    rows = [
        {
            "Compound": compound["name"],
            AFFINITY_COLUMN: round_for_display(
                (compound.get("docking") or {}).get("best_affinity_kcal_mol"), config.AFFINITY_DECIMALS
            ),
            "Docking status": _docking_status(compound),
            "Poses": (compound.get("docking") or {}).get("poses_generated"),
        }
        for compound in compounds
    ]
    if sort_by_affinity:
        # Docked compounds first, most negative affinity first; the rest keep order.
        rows.sort(key=lambda row: (row[AFFINITY_COLUMN] is None, row[AFFINITY_COLUMN] or 0.0))
    return rows


def _render_docking_setup(docking_info: dict[str, Any]) -> None:
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


def _render_docking_details(last_screening: dict[str, Any]) -> None:
    st.markdown("##### Computational Docking Results")
    docking_info = last_screening["docking"]
    if docking_info["error"]:
        st.warning(f"Docking could not run: {docking_info['error']}", icon=":material/warning:")
        return
    st.caption(AFFINITY_EXPLANATION)
    sort_choice = st.radio(
        "Sort docking results", (SORT_BY_AFFINITY, SORT_BY_DATASET), horizontal=True, label_visibility="collapsed"
    )
    st.dataframe(
        _docking_rows(last_screening["compounds"], sort_choice == SORT_BY_AFFINITY),
        column_config={
            AFFINITY_COLUMN: st.column_config.NumberColumn(AFFINITY_COLUMN, help=AFFINITY_HELP, format=AFFINITY_FORMAT),
            "Poses": st.column_config.NumberColumn("Poses", help="Poses returned by Vina.", format="%d"),
        },
        hide_index=True,
        width="stretch",
    )
    with st.expander("How to read the docking affinity"):
        st.markdown(SCORE_GUIDANCE)
    _render_docking_setup(docking_info)


def _render_property_details(compounds: list[dict[str, Any]]) -> None:
    st.markdown("##### Molecular Properties")
    st.caption(
        "Calculated by RDKit. Hover over a column heading for a short explanation. These are "
        "properties of each molecule, not evidence of activity against any virus."
    )
    st.dataframe(
        [
            {
                "Compound": compound["name"],
                **_descriptor_values(compound),
                "Status": "Calculated" if compound["valid"] else f"Failed: {compound['error']}",
            }
            for compound in compounds
        ],
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
                    ranking=compound.get("ranking"),
                )


def _render_details(last_screening: dict[str, Any]) -> None:
    _render_docking_details(last_screening)
    st.divider()
    _render_property_details(last_screening["compounds"])


# --- Report --------------------------------------------------------------------------

def _render_report(last_screening: dict[str, Any]) -> None:
    st.write(
        "A reproducible summary of this run, built directly from the pipeline results "
        "with fixed templates (no language model). It documents computational results only."
    )
    render_report(last_screening.get("report"))


# --- Page ----------------------------------------------------------------------------

RENDERERS = {
    "Overview": _render_overview,
    "Explainability": _render_explainability,
    "Ranking Method": _render_ranking_method,
    "Details": _render_details,
    "Report": _render_report,
}


def render() -> None:
    """Render the results page."""
    st.title("Results")
    last_screening = st.session_state.get(LAST_SCREENING_KEY)
    if last_screening is None:
        st.info(
            "No screening has been run yet. Open **Screening** and start a run to see results here.",
            icon=":material/hourglass_empty:",
        )
        return

    run_id = (last_screening.get("report") or {}).get("run_id")
    st.caption(f"Latest run: {last_screening['request']['target']}" + (f" · run `{run_id}`" if run_id else ""))
    for tab, title in zip(st.tabs(TAB_ORDER), TAB_ORDER):
        with tab:
            RENDERERS[title](last_screening)
