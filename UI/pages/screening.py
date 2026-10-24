"""Screening page: configure a run and send it to the backend pipeline."""

import logging

import streamlit as st

import config
from backend.screening import run_screening
from UI.utils.helpers import (
    format_status,
    is_valid_compound_count,
    render_disclaimer,
    render_pipeline,
)

logger = logging.getLogger(__name__)

# Key used to share the latest screening response with the Results page.
LAST_SCREENING_KEY = "last_screening"


def _handle_submission(target: str, n_compounds: int) -> None:
    """Validate the form input, call the backend and show the outcome."""
    if not is_valid_compound_count(n_compounds):
        st.error(
            f"Please choose between {config.MIN_COMPOUNDS} and "
            f"{config.MAX_COMPOUNDS} compounds."
        )
        return

    try:
        response = run_screening(target, n_compounds)
    except ValueError as error:
        logger.error("Invalid screening request: %s", error)
        st.error(f"Screening could not start: {error}")
        return
    except Exception:
        logger.exception("Unexpected error while starting screening")
        st.error("Something went wrong while starting the screening. Please try again.")
        return

    st.session_state[LAST_SCREENING_KEY] = response
    st.info(response["message"], icon=":material/info:")

    request = response["request"]
    target_column, count_column, status_column = st.columns(3)
    target_column.metric("Target", request["target"])
    count_column.metric("Compounds requested", request["n_compounds"])
    status_column.metric("Pipeline status", format_status(response["status"]))


def render() -> None:
    """Render the screening page."""
    st.title("Screening")
    st.write(
        "Choose a viral protein target and how many compounds to screen. "
        "Fairouz will rank the compounds by predicted research priority."
    )

    with st.form("screening_form"):
        target = st.selectbox(
            "Select Viral Target",
            config.AVAILABLE_TARGETS,
            help="Prepared protein structures the compounds will be screened against.",
        )
        n_compounds = st.slider(
            "Number of Compounds",
            min_value=config.MIN_COMPOUNDS,
            max_value=config.MAX_COMPOUNDS,
            value=config.DEFAULT_COMPOUNDS,
            help="How many compounds from the library to include in this run.",
        )
        submitted = st.form_submit_button("Start Screening", type="primary")

    if submitted:
        _handle_submission(target, n_compounds)

    st.subheader("How the pipeline will work")
    st.caption("Each stage will be implemented in a later development step.")
    render_pipeline(config.PIPELINE_STAGES)

    render_disclaimer()
