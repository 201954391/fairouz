"""Small, reusable UI helpers.

Formatting, validation and shared display blocks only. Scientific logic
belongs in ``backend/``.
"""

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
    .fairouz-step-number { font-size: 0.75rem; opacity: 0.6; letter-spacing: 0.08em; }
</style>
"""


def apply_base_styles() -> None:
    """Inject the small amount of custom CSS shared by every page."""
    st.markdown(_BASE_CSS, unsafe_allow_html=True)


def format_value(value: float | int | str | None, decimals: int = 2, unit: str = "") -> str:
    """Format a value for display, showing a dash when it is missing.

    Args:
        value: The value to format. ``None`` means "not available".
        decimals: Decimal places used for floats.
        unit: Optional unit appended after the value.

    Returns:
        A display string such as ``"-7.40 kcal/mol"`` or ``"—"``.
    """
    if value is None:
        return EMPTY_VALUE
    text = f"{value:.{decimals}f}" if isinstance(value, float) else str(value)
    return f"{text} {unit}".strip()


def format_descriptor(key: str, value: float | int | None) -> str:
    """Format a descriptor value with its unit from ``config.DESCRIPTORS``."""
    return format_value(value, unit=config.DESCRIPTORS[key]["unit"])


def descriptor_label(key: str) -> str:
    """Return a descriptor's display label, including its unit if it has one."""
    descriptor = config.DESCRIPTORS[key]
    return f"{descriptor['label']} ({descriptor['unit']})" if descriptor["unit"] else descriptor["label"]


def format_affinity(value: float | None) -> str:
    """Format a Vina affinity to one decimal place, its meaningful precision."""
    return format_value(value, decimals=1, unit="kcal/mol")


def format_status(status: str) -> str:
    """Turn a machine status such as ``"not_implemented"`` into ``"Not implemented"``."""
    readable = status.replace("_", " ").strip()
    return readable.capitalize() if readable else EMPTY_VALUE


def is_valid_compound_count(n_compounds: int) -> bool:
    """Return True if the compound count is within the configured limits."""
    return config.MIN_COMPOUNDS <= n_compounds <= config.MAX_COMPOUNDS


def render_disclaimer() -> None:
    """Display the standard research-only disclaimer."""
    st.warning(f"**Research use only.** {config.RESEARCH_DISCLAIMER}", icon=":material/science:")


def render_project_status() -> None:
    """Display the current status of each part of the project."""
    for component, status in config.PROJECT_STATUS.items():
        st.markdown(f"**{component}** · {status}")


def render_pipeline(stages: Sequence[tuple[str, str, str]]) -> None:
    """Display the screening pipeline as numbered steps with their build status.

    Args:
        stages: ``(name, description, status)`` tuples in pipeline order.
    """
    columns = st.columns(len(stages))
    for number, (column, (name, description, status)) in enumerate(zip(columns, stages), start=1):
        with column, st.container(border=True):
            st.markdown(
                f'<div class="fairouz-step-number">STEP {number} · {status.upper()}</div>',
                unsafe_allow_html=True,
            )
            st.markdown(f"**{name}**")
            st.caption(description)
