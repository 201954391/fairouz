"""Shared sidebar: branding, navigation, system status and disclaimer."""

from collections.abc import Sequence

import streamlit as st

import config
from UI.utils.helpers import render_project_status


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
        st.markdown("**System status**")
        st.caption(config.DEVELOPMENT_STAGE)
        render_project_status()

        st.divider()
        st.caption(
            ":material/science: **Research use only.** Outputs are computational "
            "predictions and require laboratory validation."
        )
        st.caption(f"v{config.PROJECT_VERSION}")
