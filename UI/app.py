"""Fairouz Streamlit entry point.

Run from the project root::

    streamlit run UI/app.py

This module only coordinates the app: page setup, navigation and the shared
sidebar. Page content lives in ``UI/pages``.
"""

import sys
from pathlib import Path

# Streamlit only puts this script's folder (UI/) on the import path. Add the
# project root so ``config``, ``backend`` and ``UI`` can be imported.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import streamlit as st  # noqa: E402

import config  # noqa: E402
from UI.components.sidebar import render_sidebar  # noqa: E402
from UI.pages import results, screening  # noqa: E402
from UI.utils.helpers import (  # noqa: E402
    apply_base_styles,
    render_disclaimer,
    render_pipeline,
    render_project_status,
)


def render_home() -> None:
    """Render the home page: what Fairouz is and how to use it."""
    st.markdown(
        f'<p class="fairouz-wordmark">{config.PROJECT_NAME.upper()}</p>'
        f'<p class="fairouz-subtitle">{config.PROJECT_SUBTITLE}</p>',
        unsafe_allow_html=True,
    )
    st.write(config.PROJECT_DESCRIPTION)

    what_column, how_column, status_column = st.columns(3)
    with what_column, st.container(border=True):
        st.markdown("**What it does**")
        st.write(
            "Screens a library of compounds against a viral protein and ranks "
            "them by predicted research priority, so lab time can focus on the "
            "most promising candidates."
        )
    with how_column, st.container(border=True):
        st.markdown("**How to use it**")
        st.write(
            "Open **Screening**, choose a target and the number of compounds, "
            "then start the run. Ranked candidates appear under **Results**."
        )
    with status_column, st.container(border=True):
        st.markdown("**Project status**")
        render_project_status()

    st.subheader("Screening pipeline")
    render_pipeline(config.PIPELINE_STAGES)

    render_disclaimer()


def main() -> None:
    """Configure the app and run the selected page."""
    config.configure_logging()
    st.set_page_config(
        page_title=f"{config.PROJECT_NAME} · {config.PROJECT_SUBTITLE}",
        page_icon="🧬",
        layout="wide",
    )
    apply_base_styles()

    # Explicit url_path values are needed because every page function is
    # called ``render``.
    pages = [
        st.Page(render_home, title="Home", icon=":material/home:", default=True),
        st.Page(screening.render, title="Screening", icon=":material/biotech:", url_path="screening"),
        st.Page(results.render, title="Results", icon=":material/analytics:", url_path="results"),
    ]
    # Hidden: our own sidebar renders the navigation links.
    current_page = st.navigation(pages, position="hidden")
    render_sidebar(pages)
    current_page.run()


if __name__ == "__main__":
    main()
