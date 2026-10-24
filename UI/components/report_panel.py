"""Report panel: status, preview and downloads for a generated report.

All report content is built by ``backend.report``; this module only shows it
and offers the files that were written.
"""

from collections.abc import Mapping
from pathlib import Path
from typing import Any

import streamlit as st

import config
from UI.utils.helpers import format_status

MIME_TYPES = {"md": "text/markdown", "json": "application/json", "html": "text/html"}
FORMAT_LABELS = {"md": "Markdown", "json": "JSON", "html": "HTML"}


def _render_status(report: Mapping[str, Any]) -> None:
    status = report["status"]
    message = f"Report status: **{format_status(status)}**"
    if report.get("run_id"):
        message += f" · Run ID `{report['run_id']}`"
    if status == config.REPORT_GENERATED:
        st.success(message, icon=":material/description:")
    else:
        st.warning(message, icon=":material/warning:")
        for problem in (report.get("errors") or {}).values():
            st.caption(problem)


def _render_preview(preview: Mapping[str, Any]) -> None:
    with st.container(border=True):
        st.markdown(f"**{preview['target']}**")
        processed, ranked, unranked = st.columns(3)
        processed.metric("Processed", preview["processed"])
        ranked.metric("Ranked", preview["ranked"])
        unranked.metric("Unranked", preview["unranked"])
        st.write(preview["outcome"])

        if preview["top"]:
            st.markdown("**Top ranked (computational screening priority)**")
            for item in preview["top"]:
                st.markdown(
                    f"- #{item['rank']} {item['name']}: Screening Priority Score {item['score']:.1f}, "
                    f"predicted Vina affinity {item['predicted_affinity_kcal_mol']:.1f} kcal/mol"
                )
        else:
            st.caption("No compound was ranked in this run, so the report has no top results.")

        st.markdown("**Key limitations**")
        for limitation in preview["key_limitations"]:
            st.markdown(f"- {limitation}")


def _render_downloads(files: Mapping[str, str]) -> None:
    if not files:
        return
    st.markdown("**Download**")
    columns = st.columns(len(files))
    for column, (fmt, path_text) in zip(columns, files.items()):
        path = Path(path_text)
        try:
            data = path.read_bytes()
        except OSError:
            column.caption(f"{FORMAT_LABELS.get(fmt, fmt)}: file no longer available")
            continue
        column.download_button(
            f"{FORMAT_LABELS.get(fmt, fmt)} report",
            data=data,
            file_name=path.name,
            mime=MIME_TYPES.get(fmt, "application/octet-stream"),
            key=f"download_report_{fmt}",
            width="stretch",
        )
    st.caption(f"Saved in `{Path(next(iter(files.values()))).parent}`.")


def render_report(report: Mapping[str, Any] | None) -> None:
    """Render the report status, preview and downloads."""
    if not report:
        st.caption("No report was produced for this run.")
        return
    _render_status(report)
    if report.get("preview"):
        _render_preview(report["preview"])
    _render_downloads(report.get("files") or {})
    if report.get("markdown"):
        with st.expander("Full report preview"):
            st.markdown(report["markdown"])
