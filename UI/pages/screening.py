"""Screening page: configure a run, show real progress, summarise the outcome."""

import logging
from typing import Any

import streamlit as st

import config
from backend.compounds import CompoundDataError
from backend.screening import run_screening
from UI.components.sidebar import vina_location
from UI.utils.helpers import (
    format_seconds,
    is_valid_compound_count,
    library_size,
    render_pipeline_status,
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
            seconds = response["timings"]["total"]
            status_box.update(label=f"Screening finished in {format_seconds(seconds)}", state="complete", expanded=False)

    if error_message:
        st.error(error_message)
        return None
    return response


def _render_timings(response: dict[str, Any]) -> None:
    timings = response["timings"]
    with st.expander("Stage timings (technical)"):
        st.dataframe(
            [
                {"Stage": stage["label"], "Time": format_seconds(timings.get(stage["key"])), "Detail": stage["detail"]}
                for stage in response["pipeline"]
            ],
            hide_index=True,
            width="stretch",
        )
        st.caption("Wall-clock times from this run on this machine; they vary between runs.")


def _render_outcome(response: dict[str, Any]) -> None:
    """Show the outcome, stage statuses and counts of a finished run."""
    if all(stage["status"] == config.STATUS_COMPLETED for stage in response["pipeline"]):
        st.success(response["message"], icon=":material/check_circle:")
    else:
        st.warning(response["message"], icon=":material/warning:")

    summary = response["summary"]
    st.markdown(
        f"**{summary['processed']} compounds processed in {format_seconds(response['timings']['total'])}.**"
    )
    render_pipeline_status(response["pipeline"])

    processed, docked, ranked, failed = st.columns(4)
    processed.metric("Processed", summary["processed"])
    docked.metric("Docked", summary["docked"], help="Compounds AutoDock Vina docked successfully.")
    ranked.metric("Ranked", summary["ranked"], help="Compounds given a Screening Priority Score.")
    failed.metric(
        "Failed",
        summary["invalid"] + summary["docking_failed"],
        help="Invalid molecules plus compounds whose ligand preparation or docking failed.",
    )

    # Problems that block whole stages are already in the message; keep other notes tidy.
    if response["warnings"]:
        with st.expander(f"Run notes ({len(response['warnings'])})"):
            for warning in response["warnings"]:
                st.markdown(f"- {warning}")
    _render_timings(response)
    st.caption("Open **Results** in the sidebar for the ranking, explanations, details and report.")


def _handle_submission(target: str, n_compounds: int) -> None:
    """Validate the form input, run the pipeline and show the outcome."""
    if not is_valid_compound_count(n_compounds):
        st.error(f"Please choose between {config.MIN_COMPOUNDS} and {config.MAX_COMPOUNDS} compounds.")
        return

    response = _run_pipeline(target, n_compounds)
    if response is None:
        return
    st.session_state[LAST_SCREENING_KEY] = response
    _render_outcome(response)


def _render_readiness(available: int | None) -> None:
    target, compounds, engine = st.columns(3)
    target.metric("Target (PDB)", config.TARGET_ID, help=config.TARGET_LABEL)
    compounds.metric("Compounds available", available if available is not None else "—")
    if vina_location():
        engine.metric("Docking engine", "Ready", help="AutoDock Vina executable found.")
    else:
        engine.metric("Docking engine", "Not found", help="Docking and ranking will be skipped.")
        st.warning(
            "AutoDock Vina was not found, so this run will calculate molecular descriptors only; "
            "docking and ranking will be skipped. See README › Installation.",
            icon=":material/warning:",
        )


def render() -> None:
    """Render the screening page."""
    st.title("Screening")
    st.write(
        "Fairouz calculates molecular descriptors with RDKit, builds a 3D model of each molecule, "
        "docks it into the target with AutoDock Vina, then ranks, explains and reports the results."
    )
    available = library_size()
    _render_readiness(available)

    with st.form("screening_form"):
        target = st.selectbox(
            "Viral target",
            config.AVAILABLE_TARGETS,
            help="The prepared protein structure the compounds are docked into.",
        )
        n_compounds = st.slider(
            "Number of compounds",
            min_value=config.MIN_COMPOUNDS,
            max_value=config.MAX_COMPOUNDS,
            value=config.DEFAULT_COMPOUNDS,
            help="If you request more than are available, all available compounds are used, without duplicates.",
        )
        submitted = st.form_submit_button("Start screening", type="primary")
    if available is not None:
        st.caption(
            f"The library holds {available} compounds; larger requests process all {available}. "
            "Docking takes a few seconds per compound."
        )

    if submitted:
        _handle_submission(target, n_compounds)

