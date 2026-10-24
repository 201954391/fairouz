"""Reusable card for displaying a single compound result."""

from collections.abc import Mapping

import streamlit as st

from UI.utils.helpers import format_value


def render_result_card(
    rank: int | None = None,
    compound_name: str | None = None,
    docking_score: float | None = None,
    ai_score: float | None = None,
    properties: Mapping[str, float | int | None] | None = None,
) -> None:
    """Display one compound result as a bordered card.

    Every argument is optional. Missing values are shown as a dash, so the
    card never displays invented numbers.

    Args:
        rank: Position of the compound in the ranked list (1 = highest priority).
        compound_name: Name or identifier of the compound.
        docking_score: Predicted binding affinity in kcal/mol (lower suggests
            stronger predicted binding).
        ai_score: Priority score from the ranking model.
        properties: Molecular properties, e.g. ``{"Molecular weight": 180.16}``.
    """
    with st.container(border=True):
        rank_label = f"#{rank} · " if rank is not None else ""
        st.markdown(f"**{rank_label}{compound_name or 'Compound not available'}**")

        docking_column, ai_column = st.columns(2)
        docking_column.metric(
            "Docking score",
            format_value(docking_score, unit="kcal/mol"),
            help="Predicted binding affinity from molecular docking. An estimate, not a measurement.",
        )
        ai_column.metric(
            "AI priority score",
            format_value(ai_score),
            help="Model-predicted research priority. Not a measure of clinical effectiveness.",
        )

        if properties:
            st.markdown("**Molecular properties**")
            property_columns = st.columns(len(properties))
            for column, (name, value) in zip(property_columns, properties.items()):
                column.metric(name, format_value(value))
        else:
            st.caption("Molecular properties: not available")

        st.caption("Computational prediction only. Requires laboratory validation.")
