"""Screening page: configure a run and send it to the backend pipeline."""

import logging
from typing import Any

import streamlit as st

import config
from backend.compounds import CompoundDataError
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


def _run_pipeline(target: str, n_compounds: int) -> dict[str, Any] | None:
    """Run the backend with live progress, returning None if it failed."""
    error_message = None
    with st.status("Running screening...", expanded=True) as status_box:
        try:
            response = run_screening(target, n_compounds, progress=status_box.write)
        except (ValueError, CompoundDataError) as error:
            logger.error("Screening could not run: %s", error)
            error_message = f"Screening could not run: {error}"
        except Exception:
            logger.exception("Unexpected error during screening")
            error_message = "Something went wrong during screening. Please try again."

        if error_message:
            status_box.update(label="Screening failed", state="error")
        else:
            status_box.update(label="Screening complete", state="complete", expanded=False)

    if error_message:
        st.error(error_message)
        return None
    return response


def _render_outcome(response: dict[str, Any]) -> None:
    """Show the counts, warnings and stage statuses of a finished run."""
    if response["docking"]["error"] is None:
        st.success(response["message"], icon=":material/check_circle:")
    else:
        st.warning(response["message"], icon=":material/warning:")

    summary = response["summary"]
    requested, loaded, analysed = st.columns(3)
    requested.metric("Requested", summary["requested"])
    loaded.metric("Loaded", summary["loaded"], help="Usable compounds in the dataset.")
    analysed.metric("Descriptors calculated", summary["valid"])
    prepared, docked, failed = st.columns(3)
    prepared.metric("3D ligands prepared", summary["ligands_prepared"])
    docked.metric("Successfully docked", summary["docked"])
    failed.metric(
        "Failed",
        summary["invalid"] + summary["docking_failed"],
        help="Invalid molecules plus compounds whose ligand preparation or docking failed.",
    )

    for warning in response["warnings"]:
        st.warning(warning, icon=":material/warning:")

    stage_lines = [
        f"- **{stage.replace('_', ' ').capitalize()}**: {format_status(status)}"
        for stage, status in response["stages"].items()
    ]
    st.markdown("\n".join(stage_lines))
    st.caption("Open **Results** in the sidebar to see the docking results and molecular properties.")


def _handle_submission(target: str, n_compounds: int) -> None:
    """Validate the form input, run the pipeline and show the outcome."""
    if not is_valid_compound_count(n_compounds):
        st.error(
            f"Please choose between {config.MIN_COMPOUNDS} and "
            f"{config.MAX_COMPOUNDS} compounds."
        )
        return

    response = _run_pipeline(target, n_compounds)
    if response is None:
        return
    st.session_state[LAST_SCREENING_KEY] = response
    _render_outcome(response)


def render() -> None:
    """Render the screening page."""
    st.title("Screening")
    st.write(
        "Choose a viral protein target and how many compounds to screen. "
        "Fairouz calculates molecular properties with RDKit, builds a 3D model "
        "of each molecule and docks it into the target with AutoDock Vina. "
        "AI ranking comes in a later stage."
    )
    st.caption(
        "Docking takes noticeably longer than property calculation: a few "
        "seconds per compound."
    )

    with st.form("screening_form"):
        target = st.selectbox(
            "Select Viral Target",
            config.AVAILABLE_TARGETS,
            help=(
                "Prepared protein structures the compounds will be docked into. "
                "6LU7 is an X-ray structure of the SARS-CoV-2 main protease."
            ),
        )
        n_compounds = st.slider(
            "Number of Compounds",
            min_value=config.MIN_COMPOUNDS,
            max_value=config.MAX_COMPOUNDS,
            value=config.DEFAULT_COMPOUNDS,
            help=(
                "How many compounds from the library to include in this run. "
                "If the library holds fewer, all of them are used."
            ),
        )
        submitted = st.form_submit_button("Start Screening", type="primary")

    if submitted:
        _handle_submission(target, n_compounds)

    st.subheader("How the pipeline works")
    render_pipeline(config.PIPELINE_STAGES)

    render_disclaimer()
