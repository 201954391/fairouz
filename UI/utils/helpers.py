"""Small, reusable UI helpers.

Formatting, validation and shared display blocks only. Scientific logic
belongs in ``backend/``.
"""

import html
import math
from collections.abc import Sequence

import streamlit as st

import config

EMPTY_VALUE = "—"
NOT_AVAILABLE_YET = "Not available yet"

AFFINITY_LABEL = "Predicted Docking Affinity"
AFFINITY_HELP = (
    "AutoDock Vina's estimate for the best predicted pose, in kcal/mol. More "
    "negative values indicate a more favourable predicted pose within this "
    "docking setup. An estimate, not a measured binding affinity."
)
PRIORITY_LABEL = "Screening Priority Score"
PRIORITY_HELP = (
    "0–100. Combines predicted docking affinity (70%) with molecular-property "
    "heuristics (30%) to prioritise compounds for further investigation. Not a "
    "probability of efficacy or success."
)
PRIORITY_EXPLANATION = (
    "Screening Priority Scores combine predicted docking affinity with "
    "molecular-property heuristics to help prioritise compounds for further "
    "investigation. They are computational screening outputs and are not "
    "evidence of efficacy, safety, or clinical benefit."
)
AFFINITY_EXPLANATION = (
    "More negative Vina scores generally indicate a more favourable predicted "
    "binding pose within this specific docking setup. These predictions are "
    "estimates and require experimental validation."
)

_BASE_CSS = """
<style>
    .block-container { padding-top: 2.5rem; max-width: 1100px; }
    .fairouz-wordmark {
        font-size: 2.6rem; font-weight: 700; letter-spacing: 0.22em;
        margin-bottom: 0;
    }
    .fairouz-subtitle { font-size: 1.15rem; opacity: 0.75; margin-top: 0; }
    .fairouz-chips { display: flex; flex-wrap: wrap; gap: 0.45rem; margin: 0.4rem 0 0.9rem; }
    .fairouz-chip {
        border: 1px solid #D9E2E7; border-radius: 0.5rem; padding: 0.35rem 0.6rem;
        font-size: 0.9rem; background: #FFFFFF; color: #1F2933; line-height: 1.35;
    }
    .fairouz-chip small { display: block; opacity: 0.7; font-size: 0.78rem; }
    .fairouz-chip-completed { border-color: #0F766E; }
    .fairouz-chip-partial { border-color: #B7791F; border-style: dashed; }
    .fairouz-chip-failed { border-color: #B42318; border-width: 2px; }
    .fairouz-chip-skipped { border-color: #98A2B3; border-style: dotted; color: #475467; }
    .fairouz-step { flex: 1 1 9rem; }
    .fairouz-step-number { font-size: 0.72rem; opacity: 0.6; letter-spacing: 0.08em; }
</style>
"""


def apply_base_styles() -> None:
    """Inject the small amount of custom CSS shared by every page."""
    st.markdown(_BASE_CSS, unsafe_allow_html=True)


def round_for_display(value: float | int | None, decimals: int) -> float | int | None:
    """Round a number for display, turning -0.0 into 0.0 (so "-0.00" never appears).

    Missing or non-finite values become None, which tables show as empty.
    """
    if value is None or isinstance(value, bool) or not math.isfinite(value):
        return None
    if decimals == 0:
        return int(round(value))
    return round(value, decimals) + 0.0


def format_value(value: float | int | str | None, decimals: int = 2, unit: str = "") -> str:
    """Format a value for display, showing a dash when it is missing.

    Args:
        value: The value to format. ``None`` means "not available".
        decimals: Decimal places used for floats.
        unit: Optional unit appended after the value.

    Returns:
        A display string such as ``"-7.4 kcal/mol"`` or ``"—"``.
    """
    if value is None:
        return EMPTY_VALUE
    if isinstance(value, float):
        rounded = round_for_display(value, decimals)
        if rounded is None:
            return EMPTY_VALUE
        text = f"{rounded:.{decimals}f}"
    else:
        text = str(value)
    return f"{text} {unit}".strip()


def format_descriptor(key: str, value: float | int | None) -> str:
    """Format a descriptor value with its configured unit and precision."""
    descriptor = config.DESCRIPTORS[key]
    return format_value(value, decimals=descriptor["decimals"], unit=descriptor["unit"])


def descriptor_number_format(key: str) -> str:
    """printf-style format for a descriptor's table column, e.g. ``"%.1f"``."""
    decimals = config.DESCRIPTORS[key]["decimals"]
    return "%d" if decimals == 0 else f"%.{decimals}f"


def descriptor_label(key: str) -> str:
    """Return a descriptor's display label, including its unit if it has one."""
    descriptor = config.DESCRIPTORS[key]
    return f"{descriptor['label']} ({descriptor['unit']})" if descriptor["unit"] else descriptor["label"]


def format_affinity(value: float | None) -> str:
    """Format a Vina affinity at its meaningful precision (one decimal)."""
    return format_value(value, decimals=config.AFFINITY_DECIMALS, unit="kcal/mol")


def format_priority(value: float | None) -> str:
    """Format a Screening Priority Score as e.g. ``"54.2 / 100"``."""
    return EMPTY_VALUE if value is None else f"{format_value(float(value), config.SCORE_DECIMALS)} / 100"


def format_seconds(seconds: float | None) -> str:
    """E.g. ``"2.4 s"``, or ``"< 0.1 s"`` for very short stages."""
    if seconds is None:
        return EMPTY_VALUE
    return "< 0.1 s" if seconds < 0.05 else f"{seconds:.1f} s"


def format_status(status: str) -> str:
    """Turn a machine status such as ``"not_generated"`` into ``"Not generated"``."""
    readable = status.replace("_", " ").strip()
    return readable.capitalize() if readable else EMPTY_VALUE


def is_valid_compound_count(n_compounds: int) -> bool:
    """Return True if the compound count is within the configured limits."""
    return config.MIN_COMPOUNDS <= n_compounds <= config.MAX_COMPOUNDS


@st.cache_data(show_spinner=False)
def library_size() -> int | None:
    """Number of usable compounds in the dataset, or None if it cannot be read.

    Cached: the dataset does not change while the app runs.
    """
    from backend.compounds import CompoundDataError, load_compounds

    try:
        return len(load_compounds().compounds)
    except CompoundDataError:
        return None


@st.cache_data(show_spinner=False)
def target_provenance() -> dict:
    """Header facts from the raw target file (method, resolution, organism), cached."""
    from backend.receptor import read_structure_provenance

    try:
        return read_structure_provenance(config.RAW_TARGET_FILE)
    except OSError:
        return {}


def render_disclaimer() -> None:
    """Display the standard research-only disclaimer."""
    st.warning(f"**Research use only.** {config.RESEARCH_DISCLAIMER}", icon=":material/science:")


def render_project_components() -> None:
    """Display the components behind each result."""
    for component, description in config.PROJECT_COMPONENTS.items():
        st.markdown(f"**{component}** · {description}")


def render_limitations() -> None:
    """Display the short limitations list."""
    for limitation in config.LIMITATIONS_SUMMARY:
        st.markdown(f"- {limitation}")


def render_workflow(stages: Sequence[tuple[str, str]]) -> None:
    """Display the workflow as numbered steps that wrap on narrow screens.

    Args:
        stages: ``(name, description)`` pairs in workflow order.
    """
    steps = "".join(
        f'<div class="fairouz-chip fairouz-step"><div class="fairouz-step-number">STEP {number}</div>'
        f"<strong>{html.escape(name)}</strong><small>{html.escape(description)}</small></div>"
        for number, (name, description) in enumerate(stages, start=1)
    )
    st.markdown(f'<div class="fairouz-chips">{steps}</div>', unsafe_allow_html=True)


# Symbol and word for each stage status. Shape and text, not colour alone,
# carry the meaning.
STAGE_STATUS_DISPLAY = {
    config.STATUS_COMPLETED: ("✓", "Completed"),
    config.STATUS_PARTIAL: ("◐", "Partial"),
    config.STATUS_FAILED: ("✕", "Failed"),
    config.STATUS_SKIPPED: ("⤼", "Skipped"),
}


def pipeline_status_html(pipeline: Sequence[dict]) -> str:
    """HTML chips for each stage: symbol, label, the status word unless completed,
    and the stage detail as a tooltip. All text is escaped."""
    chips = []
    for stage in pipeline:
        symbol, word = STAGE_STATUS_DISPLAY.get(stage["status"], ("?", format_status(stage["status"])))
        suffix = "" if stage["status"] == config.STATUS_COMPLETED else f" ({word})"
        tooltip = html.escape(f"{word}: {stage['detail']}", quote=True)
        chips.append(
            f'<span class="fairouz-chip fairouz-chip-{html.escape(stage["status"])}" title="{tooltip}">'
            f"{symbol} {html.escape(stage['label'])}{html.escape(suffix)}</span>"
        )
    return f'<div class="fairouz-chips" role="list">{"".join(chips)}</div>'


def render_pipeline_status(pipeline: Sequence[dict]) -> None:
    """Display each pipeline stage's status; hover a stage for its detail."""
    st.markdown(pipeline_status_html(pipeline), unsafe_allow_html=True)
