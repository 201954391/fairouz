"""Explainability panel: shows how one compound's ranking was calculated.

All content comes from ``backend.explainability``; this module only lays it
out. Colours are neutral (no red/green "safe/unsafe" signalling).
"""

from collections.abc import Mapping, Sequence
from typing import Any

import altair as alt
import streamlit as st

EARNED_COLOUR = "#0F766E"
NOT_EARNED_COLOUR = "#D9E2E7"
PROPERTY_COLOUR = "#64748B"

CHART_NOTE = (
    "Bars show points earned out of each component's maximum in Fairouz's formula. "
    "They are not confidence intervals or probabilities."
)
SCIENTIFIC_NOTE = (
    "This explains how Fairouz calculated the ranking. It does not describe "
    "biological efficacy, safety, toxicity or mechanism of action, and the score "
    "is not evidence of clinical benefit."
)


def _contribution_chart(rows: Sequence[Mapping[str, Any]]) -> alt.Chart:
    """Stacked horizontal bar: earned vs not-earned points per component."""
    return (
        alt.Chart(alt.Data(values=list(rows)))
        .mark_bar()
        .encode(
            x=alt.X("sum(points):Q", title="Points", scale=alt.Scale(domain=[0, 70])),
            y=alt.Y(
                "component:N",
                title=None,
                sort=["Docking", "Molecular properties"],
                axis=alt.Axis(labelOverlap=False, labelLimit=180),
            ),
            color=alt.Color(
                "segment:N",
                title=None,
                scale=alt.Scale(domain=["Earned", "Not earned"], range=[EARNED_COLOUR, NOT_EARNED_COLOUR]),
                legend=alt.Legend(orient="bottom"),
            ),
            order=alt.Order("segment:N", sort="ascending"),
            tooltip=[
                alt.Tooltip("component:N", title="Component"),
                alt.Tooltip("segment:N", title="Points"),
                alt.Tooltip("points:Q", title="Value", format=".1f"),
            ],
        )
        .properties(height=alt.Step(34))  # One fixed-height row per component.
    )


def overview_chart(explanations: Sequence[Mapping[str, Any]]) -> alt.Chart:
    """Stacked bar per ranked compound: docking points plus property points."""
    rows = []
    for explanation in explanations:
        if not explanation["ranked"]:
            continue
        label = f"#{explanation['rank']} {explanation['name']}"
        rows.append({"compound": label, "rank": explanation["rank"], "part": "Docking", "points": explanation["docking"]["points"]})
        rows.append(
            {"compound": label, "rank": explanation["rank"], "part": "Molecular properties", "points": explanation["properties"]["points"]}
        )
    return (
        alt.Chart(alt.Data(values=rows))
        .mark_bar()
        .encode(
            x=alt.X("sum(points):Q", title="Screening Priority Score (points)", scale=alt.Scale(domain=[0, 100])),
            y=alt.Y(
                "compound:N",
                title=None,
                sort=alt.EncodingSortField(field="rank", order="ascending"),
                axis=alt.Axis(labelOverlap=False, labelLimit=180),
            ),
            color=alt.Color(
                "part:N",
                title=None,
                scale=alt.Scale(domain=["Docking", "Molecular properties"], range=[EARNED_COLOUR, PROPERTY_COLOUR]),
                legend=alt.Legend(orient="bottom"),
            ),
            order=alt.Order("part:N", sort="ascending"),
            tooltip=[
                alt.Tooltip("compound:N", title="Compound"),
                alt.Tooltip("part:N", title="Component"),
                alt.Tooltip("points:Q", title="Points", format=".1f"),
            ],
        )
        .properties(height=alt.Step(26))  # A fixed row height, so every compound gets a label.
    )


def _descriptor_rows(items: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "Property": item["label"],
            "Value": item["display_value"],
            "Configured range": item["range_text"],
            "Position": item["position"].capitalize(),
            "Sub-score": item["subscore"],
            "Points": item["points"],
            "Interpretation": item["interpretation"],
        }
        for item in items
    ]


def _render_summary(explanation: Mapping[str, Any]) -> None:
    docking, properties = explanation["docking"], explanation["properties"]
    with st.container(border=True):
        st.markdown(f"#### #{explanation['rank']} · {explanation['name']}")
        left, right = st.columns(2)
        left.markdown(f"Screening Priority Score: **{explanation['score_display']} / 100**")
        left.markdown(f"Primary driver: **{explanation['primary_driver']}**", help=explanation["driver_note"])
        right.markdown(f"Docking contribution: **{docking['points']:.1f} / {docking['max_points']:.0f}**")
        right.markdown(f"Property contribution: **{properties['points']:.1f} / {properties['max_points']:.0f}**")
        st.write(explanation["summary"])


def _render_ranked(explanation: Mapping[str, Any]) -> None:
    _render_summary(explanation)

    st.markdown("**Contribution breakdown**")
    st.altair_chart(_contribution_chart(explanation["chart_rows"]), width="stretch")
    st.caption(CHART_NOTE)
    st.write(explanation["docking"]["text"])
    st.caption(explanation["docking"]["note"])

    st.markdown("**Descriptor breakdown**")
    st.dataframe(
        _descriptor_rows(explanation["properties"]["items"]),
        column_config={
            "Sub-score": st.column_config.ProgressColumn("Sub-score", min_value=0, max_value=1, format="%.2f"),
            "Points": st.column_config.NumberColumn(
                "Points", help=f"Out of {explanation['properties']['items'][0]['max_points']:.0f} per property.", format="%.1f"
            ),
        },
        hide_index=True,
        width="stretch",
    )

    strengths_column, limitations_column = st.columns(2)
    with strengths_column:
        st.markdown("**Main strengths**")
        st.markdown(f"*Main positive factor:* {explanation['main_positive_factor']}")
        for strength in explanation["other_strengths"]:
            st.markdown(f"- {strength}")
    with limitations_column:
        st.markdown("**Main limitations**")
        st.markdown(f"*Main limiting factor:* {explanation['main_limiting_factor']}.")
        for limitation in explanation["other_limitations"]:
            st.markdown(f"- {limitation}")

    st.markdown("**Why this rank?**")
    for comparison in (explanation["comparison_above"], explanation["comparison_below"]):
        if comparison:
            st.markdown(f"- {comparison}")
    for note in explanation["tie_notes"]:
        st.info(f"Tie-break: {note}", icon=":material/balance:")


def _render_unranked(explanation: Mapping[str, Any]) -> None:
    with st.container(border=True):
        st.markdown(f"#### Not ranked · {explanation['name']}")
        st.markdown(f"Status: **{explanation['status']}** · Pipeline stage: **{explanation['failed_stage']}**")
        st.write(explanation["summary"])
        with st.expander("Details"):
            st.write(explanation["detail"])


def render_explanation(explanation: Mapping[str, Any]) -> None:
    """Render the explanation for one compound, ranked or not."""
    if explanation["ranked"]:
        _render_ranked(explanation)
    else:
        _render_unranked(explanation)
    st.caption(SCIENTIFIC_NOTE)
