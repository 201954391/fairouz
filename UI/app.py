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
from UI.pages import home, results, screening  # noqa: E402
from UI.utils.helpers import apply_base_styles  # noqa: E402


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
        st.Page(home.render, title="Home", icon=":material/home:", default=True),
        st.Page(screening.render, title="Screening", icon=":material/biotech:", url_path="screening"),
        st.Page(results.render, title="Results", icon=":material/analytics:", url_path="results"),
    ]
    home.screening_page = pages[1]
    # Hidden: our own sidebar renders the navigation links.
    current_page = st.navigation(pages, position="hidden")
    render_sidebar(pages)
    current_page.run()


if __name__ == "__main__":
    main()
