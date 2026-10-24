"""Shared sidebar: branding, navigation, components and docking readiness."""

from collections.abc import Sequence

import streamlit as st

import config
from backend.docking import VinaNotFoundError, find_vina_executable
from UI.utils.helpers import render_project_components


@st.cache_data(ttl=60, show_spinner=False)
def vina_location() -> str | None:
    """Where AutoDock Vina was found, or None. Re-checked at most once a minute."""
    try:
        return str(find_vina_executable())
    except VinaNotFoundError:
        return None


def _render_docking_readiness() -> None:
    location = vina_location()
    if location:
        st.caption(":material/check_circle: **AutoDock Vina found**", help=location)
    else:
        st.caption(
            ":material/cancel: **AutoDock Vina not found**: docking and ranking will be skipped. "
            "See README › Installation."
        )


def render_sidebar(pages: Sequence[st.Page]) -> None:
    """Render the Fairouz sidebar.

    Args:
        pages: The app's pages, in the order they should appear.
    """
    with st.sidebar:
        st.markdown(f"### {config.PROJECT_NAME.upper()}")
        st.caption(config.PROJECT_SUBTITLE)

        st.divider()
        for page in pages:
            st.page_link(page)

        st.divider()
        render_project_components()
        _render_docking_readiness()

        st.divider()
        st.caption(
            ":material/science: **Research use only.** Outputs are computational "
            "and require experimental validation."
        )
        st.caption(f"v{config.PROJECT_VERSION} · {config.DEVELOPMENT_STAGE}")
