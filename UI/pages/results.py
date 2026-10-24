"""Results page: where screening outputs will be displayed."""

import streamlit as st

from UI.components.result_card import render_result_card
from UI.pages.screening import LAST_SCREENING_KEY
from UI.utils.helpers import format_status, render_disclaimer

# (section title, what it will show) in display order.
RESULT_SECTIONS = (
    (
        "Top Ranked Compounds",
        "The compounds Fairouz predicts are most worth investigating further, "
        "with the main reasons each one ranked highly.",
    ),
    (
        "Docking Results",
        "Predicted binding affinities (kcal/mol) from molecular docking. "
        "These are estimates of binding, not measurements.",
    ),
    (
        "Molecular Properties",
        "Calculated properties such as molecular weight, LogP, TPSA and "
        "hydrogen bond donors/acceptors.",
    ),
    (
        "AI Ranking",
        "How the ranking model combined these features, including the "
        "uncertainty of each prediction.",
    ),
    (
        "Final Report",
        "A downloadable research summary covering the target, data sources, "
        "top candidates, limitations and next steps for laboratory validation.",
    ),
)


def _render_last_request() -> None:
    """Show the most recent screening request from this session, if any."""
    last_screening = st.session_state.get(LAST_SCREENING_KEY)
    if last_screening is None:
        st.info(
            "No screening has been run yet. Results will appear here after the "
            "screening pipeline is implemented.",
            icon=":material/hourglass_empty:",
        )
        return

    request = last_screening["request"]
    st.info(
        f"Last request: **{request['target']}**, {request['n_compounds']} compounds. "
        f"Status: **{format_status(last_screening['status'])}**. "
        "No results were generated.",
        icon=":material/hourglass_empty:",
    )


def render() -> None:
    """Render the results page."""
    st.title("Results")
    st.write("Ranked compounds, docking results and the final report will appear here.")
    _render_last_request()

    tabs = st.tabs([title for title, _ in RESULT_SECTIONS])
    for tab, (title, description) in zip(tabs, RESULT_SECTIONS):
        with tab:
            st.markdown(f"#### {title}")
            st.write(description)
            st.caption("No data available yet.")

    # Show the card layout with empty values so the design can be reviewed.
    with tabs[0]:
        st.markdown("**Result card layout preview**")
        render_result_card()

    render_disclaimer()
