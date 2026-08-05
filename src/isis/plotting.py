"""
Publication-style plotting for ISIS predictions.

A single place that defines the visual language for every ISIS figure, so
prediction profiles and benchmark summaries share one set of axes, fonts,
colours and spacing rather than each script inventing its own.

Design rules enforced here (deliberately, not by taste):

- One y-axis per panel. Measures on different scales become separate stacked
  panels sharing the residue axis, never a twinned second axis.
- Categorical hues are assigned from a fixed, validated order and never
  cycled or generated. The order passes colour-vision-deficiency separation
  and normal-vision separation gates for adjacent pairs.
- Magnitude is encoded once. A per-residue score is drawn as height, in a
  single hue - never as a rainbow ramp keyed to the value, which would
  double-encode what the height already says.
- Continuous magnitude (heatmaps) uses one hue, light to dark.
- Chrome recedes: hairline horizontal gridlines, no top/right spines, muted
  tick ink, so the data is the darkest thing in the frame.
- Three of the categorical hues fall below 3:1 contrast on a light surface,
  so charts using them carry visible value labels (the "relief" rule).

Typical use:

    from isis.plotting import plot_profile, save_figure
    fig = plot_profile(sequence, results)
    save_figure(fig, "profile.png")
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# ---------------------------------------------------------------------------
# Palette
# ---------------------------------------------------------------------------

# Fixed categorical order. Assign by slot index; never reorder per-figure and
# never generate a 9th hue - fold extra series into "other" or facet instead.
SERIES: Tuple[str, ...] = (
    "#2a78d6",  # 1 blue
    "#eb6834",  # 2 orange
    "#1baf7a",  # 3 aqua
    "#eda100",  # 4 yellow
    "#e87ba4",  # 5 magenta
    "#008300",  # 6 green
    "#4a3aa7",  # 7 violet
    "#e34948",  # 8 red
)

# Hues below 3:1 against the light surface. Charts that use these must show
# visible value labels so colour is not the only channel carrying identity.
LOW_CONTRAST_SERIES = frozenset({"#1baf7a", "#eda100", "#e87ba4"})

# Single-hue sequential ramp for continuous magnitude (heatmaps).
SEQUENTIAL = ("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")

# Discrete ordinal steps, spaced for visible lightness gaps.
ORDINAL = ("#86b6ef", "#3987e5", "#256abf", "#104281")

# Diverging pair: warm/cool poles with a neutral midpoint that reads as "nothing".
DIVERGING = ("#2a78d6", "#f0efec", "#e34948")

INK = {
    "surface": "#fcfcfb",
    "primary": "#0b0b0b",
    "secondary": "#52514e",
    "muted": "#898781",
    "grid": "#e1e0d9",
    "axis": "#c3c2b7",
}

# Status colours are reserved for state, never reused as a series hue.
STATUS = {
    "good": "#0ca30c",
    "warning": "#fab219",
    "serious": "#ec835a",
    "critical": "#d03b3b",
}

SEQUENTIAL_CMAP = LinearSegmentedColormap.from_list("isis_seq", SEQUENTIAL)


def _stable_series_colors(keys: Sequence[str]) -> Dict[str, str]:
    """
    Map names to palette slots by sorted position.

    Sorting rather than call-order means a method keeps its colour whether or
    not its neighbours are present, so dropping a series from a figure never
    repaints the survivors.
    """
    return {k: SERIES[i % len(SERIES)] for i, k in enumerate(sorted(keys))}


# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------

def apply_style() -> None:
    """Install the ISIS matplotlib style. Idempotent; safe to call repeatedly."""
    plt.rcParams.update({
        "figure.facecolor": INK["surface"],
        "figure.dpi": 110,
        "savefig.dpi": 200,
        "savefig.facecolor": INK["surface"],
        "savefig.bbox": "tight",

        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial", "sans-serif"],
        "font.size": 9,

        "axes.facecolor": INK["surface"],
        "axes.edgecolor": INK["axis"],
        "axes.linewidth": 0.8,
        "axes.labelcolor": INK["secondary"],
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "axes.titlecolor": INK["primary"],
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.titlepad": 8,
        # Chrome recedes: only a hairline horizontal grid, always behind the data.
        "axes.grid": True,
        "axes.grid.axis": "y",
        "axes.axisbelow": True,
        "axes.spines.top": False,
        "axes.spines.right": False,

        "grid.color": INK["grid"],
        "grid.linewidth": 0.6,

        "xtick.color": INK["muted"],
        "ytick.color": INK["muted"],
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.direction": "out",
        "ytick.direction": "out",

        "legend.frameon": False,
        "legend.fontsize": 8,
        "legend.labelcolor": INK["secondary"],

        "lines.linewidth": 1.6,
        "lines.solid_capstyle": "round",
    })


def save_figure(fig, path: str) -> str:
    """Write a figure and close it, returning the path."""
    fig.savefig(path)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def align_to_sequence(scores: Sequence[float], seq_len: int) -> np.ndarray:
    """
    Pad a windowed score array to one value per residue.

    Sliding-window methods return fewer scores than residues (a window of w
    loses w-1 positions), and those scores describe window *centres*. Padding
    symmetrically with NaN keeps score i aligned with residue i and leaves the
    unscored ends as genuine gaps rather than fabricated zeros.
    """
    arr = np.asarray(scores, dtype=float)
    if arr.size == seq_len:
        return arr
    if arr.size > seq_len:
        return arr[:seq_len]
    total = seq_len - arr.size
    left = total // 2
    return np.concatenate([np.full(left, np.nan), arr, np.full(total - left, np.nan)])


def _spans_above(values: np.ndarray, threshold: float, min_len: int = 1
                 ) -> List[Tuple[int, int]]:
    """Contiguous 1-indexed [start, end] runs at or above a threshold."""
    spans: List[Tuple[int, int]] = []
    start: Optional[int] = None
    for i, v in enumerate(values):
        hot = not np.isnan(v) and v >= threshold
        if hot and start is None:
            start = i
        elif not hot and start is not None:
            if i - start >= min_len:
                spans.append((start + 1, i))
            start = None
    if start is not None and len(values) - start >= min_len:
        spans.append((start + 1, len(values)))
    return spans


def _draw_highlight(ax, spans: Optional[Sequence[Tuple[int, int]]]) -> None:
    """
    Shade reference spans behind the data.

    Kept semi-transparent and at the bottom of the z-order: on antigens whose
    reference epitopes cover most of the sequence, an opaque band would
    dominate the panel and bury the trace it is meant to contextualise.
    """
    if not spans:
        return
    for start, end in spans:
        ax.axvspan(start, end, color=INK["grid"], alpha=0.6, lw=0, zorder=0)


def _residue_axis(ax, seq_len: int, label: bool = False) -> None:
    ax.set_xlim(1, seq_len)
    ax.margins(x=0)
    if label:
        ax.set_xlabel("Residue position")


def _label_bars(ax, bars, values, fmt: str = "{:.2f}", fontsize: int = 7) -> None:
    """
    Print each bar's value above it.

    Required, not decorative: several palette hues sit below 3:1 contrast on a
    light surface, so the figure must not rely on colour alone to be readable.
    """
    for bar, val in zip(bars, values):
        ax.annotate(fmt.format(val),
                    xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
                    xytext=(0, 2), textcoords="offset points",
                    ha="center", va="bottom",
                    fontsize=fontsize, color=INK["secondary"])


# ---------------------------------------------------------------------------
# Prediction figures
# ---------------------------------------------------------------------------

def plot_profile(
    sequence: str,
    results: Dict[str, Any],
    title: str = "B-cell epitope prediction profile",
    subtitle: Optional[str] = None,
    highlight: Optional[Sequence[Tuple[int, int]]] = None,
) -> plt.Figure:
    """
    Per-residue score profile, one stacked panel per method.

    Small multiples rather than overlaid lines: the methods use different score
    scales (surface-accessibility ratios, hydrophilicity, propensity indices),
    and stacking them on a shared residue axis keeps each on its own honest
    y-scale while still lining positions up vertically for comparison.

    Args:
        sequence: the scored sequence, used only for its length.
        results: method name -> object with .scores and .threshold (an ISIS
            Prediction), or a dict with those keys.
        highlight: optional spans (1-indexed, inclusive) drawn behind every
            panel - e.g. experimentally confirmed epitopes.
    """
    apply_style()
    seq_len = len(sequence)
    names = sorted(results)
    colors = _stable_series_colors(names)

    fig, axes = plt.subplots(
        len(names), 1, figsize=(11, 1.35 * len(names) + 1.1),
        sharex=True, squeeze=False,
    )
    axes = axes[:, 0]

    for ax, name in zip(axes, names):
        res = results[name]
        scores = getattr(res, "scores", None)
        threshold = getattr(res, "threshold", None)
        if scores is None and isinstance(res, dict):
            scores = res.get("scores")
            threshold = res.get("threshold")

        values = align_to_sequence(scores, seq_len)
        x = np.arange(1, seq_len + 1)
        color = colors[name]

        _draw_highlight(ax, highlight)

        ax.plot(x, values, color=color, lw=1.3, zorder=3)

        n_regions = 0
        if threshold is not None:
            # Fill only where the curve exceeds its threshold, measured from
            # the threshold itself. Shading whole regions with axvspan is
            # unreadable once a method calls 70+ regions, and filling down to
            # zero would force the y-axis to include zero - flattening methods
            # whose scores sit near 1.0 and hiding the variation that matters.
            ax.fill_between(x, values, threshold,
                            where=(values >= threshold), interpolate=True,
                            color=color, alpha=0.32, lw=0, zorder=2)
            ax.axhline(threshold, color=INK["muted"], lw=0.9, ls=(0, (4, 3)),
                       zorder=4)
            n_regions = len(_spans_above(values, threshold, min_len=1))

        # Y-range from the data, with a little headroom. Never anchored to 0.
        finite = values[np.isfinite(values)]
        if finite.size:
            lo, hi = float(finite.min()), float(finite.max())
            pad = max((hi - lo) * 0.12, 1e-6)
            ax.set_ylim(lo - pad, hi + pad)
        ax.yaxis.set_major_locator(plt.MaxNLocator(3))

        # Region count rides in the axis label, where nothing can collide
        # with the trace.
        label = f"{name}\n{n_regions} regions" if threshold is not None else name
        ax.set_ylabel(label, rotation=0, ha="right", va="center",
                      labelpad=12, fontsize=8.5, color=INK["secondary"],
                      linespacing=1.5)
        _residue_axis(ax, seq_len)

    _residue_axis(axes[-1], seq_len, label=True)

    fig.suptitle(title, x=0.012, y=0.995, ha="left", fontsize=12,
                 fontweight="bold", color=INK["primary"])
    caption = subtitle or f"{seq_len} residues"
    caption += "   ·   dashed line: method threshold; fill marks residues above it"
    if highlight:
        caption += "   ·   grey bands: reference epitopes"
    fig.text(0.012, 0.958, caption, ha="left", fontsize=8.5,
             color=INK["muted"])
    fig.tight_layout(rect=(0, 0, 1, 0.945))
    return fig


def plot_consensus(
    sequence: str,
    votes: Sequence[float],
    min_methods: int = 3,
    n_methods: Optional[int] = None,
    title: str = "B-cell epitope consensus",
    highlight: Optional[Sequence[Tuple[int, int]]] = None,
) -> plt.Figure:
    """
    Agreement across methods, as a count per residue.

    Drawn as single-hue bars: the count is already carried by height, so
    colouring bars by their value would encode the same variable twice and
    burn the one free channel. Residues meeting the agreement threshold are
    picked out by a marker strip beneath the axis instead.
    """
    apply_style()
    seq_len = len(sequence)
    values = align_to_sequence(votes, seq_len)
    x = np.arange(1, seq_len + 1)

    fig, ax = plt.subplots(figsize=(11, 2.9))

    _draw_highlight(ax, highlight)

    ax.bar(x, np.nan_to_num(values, nan=0.0), width=1.0,
           color=SERIES[0], lw=0, zorder=2)
    ax.axhline(min_methods, color=INK["muted"], lw=0.9, ls=(0, (4, 3)), zorder=3)

    top = n_methods if n_methods else int(np.nanmax(values)) + 1
    top = max(top, min_methods + 1)
    ax.set_ylim(0, top)

    # Marker strip below the axis, dropped clear of the tick labels (which get
    # extra pad) so the two never interleave.
    spans = _spans_above(values, min_methods, min_len=1)
    strip_y = -top * 0.13
    for start, end in spans:
        ax.plot([start, end], [strip_y, strip_y], color=STATUS["critical"],
                lw=3.5, solid_capstyle="butt", clip_on=False, zorder=4)
    ax.tick_params(axis="x", pad=20)

    ax.set_ylabel("Methods in agreement")
    _residue_axis(ax, seq_len, label=True)

    handles = [
        Patch(facecolor=SERIES[0], label="Methods calling this residue"),
        Line2D([0], [0], color=INK["muted"], lw=0.9, ls=(0, (4, 3)),
               label=f"Consensus threshold ({min_methods})"),
        Line2D([0], [0], color=STATUS["critical"], lw=3,
               label=f"Consensus regions (n={len(spans)})"),
    ]
    if highlight:
        handles.append(Patch(facecolor=INK["grid"], label="Reference epitopes"))
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, 1.22), ncol=4)

    fig.suptitle(title, x=0.012, y=1.02, ha="left", fontsize=12,
                 fontweight="bold", color=INK["primary"])
    fig.tight_layout(rect=(0, 0, 1, 0.99))
    return fig


def plot_mhc_binding(
    sequence: str,
    peptides: Sequence[Dict[str, Any]],
    allele: str,
    strong_nm: float = 50.0,
    weak_nm: float = 500.0,
    title: Optional[str] = None,
) -> plt.Figure:
    """
    Predicted MHC binding affinity across a sequence.

    IC50 is plotted on a log axis and inverted, so stronger binders sit higher
    and the visual "up is more immunogenic" reading matches the biology.
    Affinity class is a status, not an identity, so it uses status colours.
    """
    apply_style()
    seq_len = len(sequence)

    fig, ax = plt.subplots(figsize=(11, 3.1))

    starts = np.array([p["start"] for p in peptides], dtype=float)
    ic50s = np.array([max(float(p["ic50"]), 1e-3) for p in peptides])

    ax.axhspan(strong_nm, weak_nm, color=INK["grid"], lw=0, zorder=0)

    ax.scatter(starts, ic50s, s=13, color=SERIES[0], alpha=0.55,
               lw=0, zorder=2, label="Predicted peptide")
    for cutoff, key, label in ((strong_nm, "critical", f"Strong (<{strong_nm:g} nM)"),
                               (weak_nm, "warning", f"Weak (<{weak_nm:g} nM)")):
        ax.axhline(cutoff, color=STATUS[key], lw=1.1, ls=(0, (4, 3)),
                   zorder=3, label=label)

    ax.set_yscale("log")
    ax.invert_yaxis()
    ax.set_ylabel("Predicted IC$_{50}$ (nM, log)")
    _residue_axis(ax, seq_len, label=True)

    n_strong = int((ic50s < strong_nm).sum())
    n_weak = int(((ic50s >= strong_nm) & (ic50s < weak_nm)).sum())

    ax.legend(loc="upper left", bbox_to_anchor=(0, 1.2), ncol=3)
    fig.suptitle(title or f"MHC binding prediction — {allele}",
                 x=0.012, y=1.03, ha="left", fontsize=12,
                 fontweight="bold", color=INK["primary"])
    fig.text(0.012, 0.94,
             f"{len(peptides)} peptides scored   ·   {n_strong} strong, "
             f"{n_weak} weak binders   ·   lower IC$_{{50}}$ = tighter binding (plotted upward)",
             ha="left", fontsize=8.5, color=INK["muted"])
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return fig


def plot_call_matrix(
    sequence: str,
    results: Dict[str, Any],
    title: str = "Which methods call which residues",
    highlight: Optional[Sequence[Tuple[int, int]]] = None,
) -> plt.Figure:
    """
    Methods x residues matrix of binary epitope calls.

    Deliberately binary rather than a continuous heatmap of scores. At a few
    hundred residues each cell is only a pixel or two wide, and a continuous
    map of per-residue scores at that density reads as static - the eye cannot
    recover anything from it. The called/not-called question survives the
    density, and answers what the consensus track cannot: *which* methods
    agree on a given region, not merely how many.

    Each method keeps its own categorical hue, so a row is identifiable by
    position and colour together.
    """
    apply_style()
    seq_len = len(sequence)
    names = sorted(results)
    colors = _stable_series_colors(names)
    x = np.arange(1, seq_len + 1)

    fig, ax = plt.subplots(figsize=(11, 0.36 * len(names) + 2.1))

    _draw_highlight(ax, highlight)

    for row, name in enumerate(names):
        res = results[name]
        epitopes = getattr(res, "epitopes", None)
        called = np.zeros(seq_len, dtype=bool)
        if epitopes is not None:
            for ep in epitopes:
                called[ep.start - 1:ep.end] = True
        else:
            scores = align_to_sequence(res.get("scores", []), seq_len)
            thr = res.get("threshold")
            if thr is not None:
                called = np.nan_to_num(scores, nan=-np.inf) >= thr

        # A thin band per method: bar height 0.7 leaves a surface gap between
        # rows so adjacent fills never touch.
        ax.fill_between(x, row - 0.35, row + 0.35, where=called,
                        color=colors[name], lw=0, zorder=2, step="mid")

        pct = 100.0 * called.sum() / seq_len
        ax.annotate(f"{pct:.0f}%", xy=(1.005, row), xycoords=("axes fraction", "data"),
                    va="center", ha="left", fontsize=7.5, color=INK["muted"])

    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=8.5, color=INK["secondary"])
    ax.set_ylim(len(names) - 0.5, -0.5)
    ax.grid(False)
    _residue_axis(ax, seq_len, label=True)

    fig.suptitle(title, x=0.012, y=1.005, ha="left", fontsize=12,
                 fontweight="bold", color=INK["primary"])
    caption = "Filled = residue called an epitope   ·   right margin: % of sequence called"
    if highlight:
        caption += "   ·   grey bands: reference epitopes"
    fig.text(0.012, 0.93, caption, ha="left", fontsize=8.5, color=INK["muted"])
    fig.tight_layout(rect=(0, 0, 1, 0.915))
    return fig


# ---------------------------------------------------------------------------
# Benchmark figures
# ---------------------------------------------------------------------------

METRIC_LABELS = {
    "sens": "Sensitivity",
    "spec": "Specificity",
    "ppv": "PPV",
    "f1": "F1",
    "mcc": "MCC",
}


def plot_metric_comparison(
    per_method: Dict[str, Dict[str, float]],
    metrics: Sequence[str] = ("sens", "spec", "ppv", "f1"),
    title: str = "Per-method benchmark performance",
    subtitle: Optional[str] = None,
) -> plt.Figure:
    """
    Grouped bars: one group per method, one bar per metric.

    Metrics are nominal identities, so each takes a fixed categorical slot -
    not a value ramp. Every bar is labelled with its value, which the palette's
    low-contrast hues require and which a reader wants anyway.
    """
    apply_style()
    names = sorted(per_method)
    metric_colors = {m: SERIES[i % len(SERIES)] for i, m in enumerate(metrics)}

    x = np.arange(len(names))
    width = 0.8 / len(metrics)

    fig, ax = plt.subplots(figsize=(1.6 * len(names) + 3.4, 4.2))

    for i, metric in enumerate(metrics):
        vals = [per_method[n].get(metric, np.nan) for n in names]
        offset = (i - (len(metrics) - 1) / 2) * width
        bars = ax.bar(x + offset, vals, width * 0.92,
                      color=metric_colors[metric], lw=0,
                      label=METRIC_LABELS.get(metric, metric))
        _label_bars(ax, bars, vals)

    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=20, ha="right", fontsize=8.5,
                       color=INK["secondary"])
    ax.set_ylabel("Score")
    ax.set_ylim(0, 1.0)
    ax.legend(loc="upper left", bbox_to_anchor=(0, 1.14),
              ncol=len(metrics))

    fig.suptitle(title, x=0.012, y=1.015, ha="left", fontsize=12,
                 fontweight="bold", color=INK["primary"])
    if subtitle:
        fig.text(0.012, 0.955, subtitle, ha="left", fontsize=8.5,
                 color=INK["muted"])
    fig.tight_layout(rect=(0, 0, 1, 0.945 if subtitle else 0.975))
    return fig


def plot_mcc_ranking(
    mcc_by_method: Dict[str, float],
    title: str = "Correlation with ground truth (MCC)",
    subtitle: Optional[str] = None,
) -> plt.Figure:
    """
    Horizontal MCC bars against a zero reference.

    MCC has a meaningful zero - it is exactly the score of random guessing -
    so the figure is built around that line. Bars are coloured by sign as a
    status (better/worse than chance), which is polarity rather than identity,
    with the value printed on every bar.
    """
    apply_style()
    order = sorted(mcc_by_method, key=lambda k: mcc_by_method[k])
    vals = [mcc_by_method[k] for k in order]
    colors = [STATUS["good"] if v > 0 else STATUS["critical"] for v in vals]

    fig, ax = plt.subplots(figsize=(7.6, 0.44 * len(order) + 2.0))

    y = np.arange(len(order))
    ax.barh(y, vals, height=0.62, color=colors, lw=0, zorder=2)
    ax.axvline(0, color=INK["secondary"], lw=1.0, zorder=3)

    for yi, v in zip(y, vals):
        ax.annotate(f"{v:+.3f}",
                    xy=(v, yi),
                    xytext=(4 if v >= 0 else -4, 0), textcoords="offset points",
                    ha="left" if v >= 0 else "right", va="center",
                    fontsize=7.5, color=INK["secondary"])

    ax.set_yticks(y)
    ax.set_yticklabels(order, fontsize=8.5, color=INK["secondary"])
    ax.set_xlabel("Matthews correlation coefficient")
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)

    pad = max(0.05, max(abs(min(vals)), abs(max(vals))) * 0.35)
    ax.set_xlim(min(min(vals), 0) - pad, max(max(vals), 0) + pad)

    fig.suptitle(title, x=0.012, y=1.02, ha="left", fontsize=12,
                 fontweight="bold", color=INK["primary"])
    if subtitle:
        fig.text(0.012, 0.95, subtitle, ha="left", fontsize=8.5,
                 color=INK["muted"])
    fig.tight_layout(rect=(0, 0, 1, 0.94 if subtitle else 0.98))
    return fig


def plot_threshold_sweep(
    thresholds: Sequence[float],
    metrics_by_threshold: Dict[str, Sequence[float]],
    title: str = "Effect of consensus threshold",
    subtitle: Optional[str] = None,
    xlabel: str = "Consensus threshold (methods in agreement)",
) -> plt.Figure:
    """
    Metric trajectories across a threshold sweep, all on one 0-1 axis.

    Every series here is a rate in [0, 1] except MCC, which shares the range
    but has a meaningful zero - so a zero reference line is drawn rather than
    giving MCC a second axis. Series are directly labelled at their right-hand
    end in addition to the legend, so identity never rests on colour alone.
    """
    apply_style()
    names = list(metrics_by_threshold)
    colors = {n: SERIES[i % len(SERIES)] for i, n in enumerate(names)}
    x = np.asarray(thresholds, dtype=float)

    fig, ax = plt.subplots(figsize=(8.0, 4.4))

    ax.axhline(0, color=INK["axis"], lw=0.9, zorder=1)

    ends = []
    for name in names:
        vals = np.asarray(metrics_by_threshold[name], dtype=float)
        ax.plot(x, vals, color=colors[name], lw=1.8, marker="o",
                markersize=4.5, markeredgecolor=INK["surface"],
                markeredgewidth=0.8, zorder=3,
                label=METRIC_LABELS.get(name, name))
        ends.append([float(vals[-1]), name])

    # Nudge end labels apart where series converge: several metrics land near
    # the same value at high thresholds, and stacked text is unreadable.
    ends.sort()
    span = max(e[0] for e in ends) - min(e[0] for e in ends)
    min_gap = max(span * 0.055, 0.035)
    for i in range(1, len(ends)):
        if ends[i][0] - ends[i - 1][0] < min_gap:
            ends[i][0] = ends[i - 1][0] + min_gap
    for y_pos, name in ends:
        ax.annotate(METRIC_LABELS.get(name, name),
                    xy=(x[-1], y_pos), xytext=(7, 0),
                    textcoords="offset points", va="center",
                    fontsize=7.5, color=colors[name])

    ax.set_xticks(x)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Score")
    ax.set_xlim(x.min() - 0.15, x.max() + (x.max() - x.min()) * 0.22)
    ax.legend(loc="upper left", bbox_to_anchor=(0, 1.13), ncol=len(names))

    fig.suptitle(title, x=0.012, y=1.015, ha="left", fontsize=12,
                 fontweight="bold", color=INK["primary"])
    if subtitle:
        fig.text(0.012, 0.95, subtitle, ha="left", fontsize=8.5,
                 color=INK["muted"])
    fig.tight_layout(rect=(0, 0, 1, 0.94 if subtitle else 0.975))
    return fig


def plot_score_distribution(
    values: Sequence[float],
    xlabel: str = "F1 score",
    title: str = "Distribution across antigens",
    subtitle: Optional[str] = None,
    bins: int = 16,
) -> plt.Figure:
    """
    Histogram of a per-antigen metric, with mean and median marked.

    A single distribution is one series, so it takes slot 1 and needs no
    legend box for the bars - the title names them. The two reference lines
    are labelled directly.
    """
    apply_style()
    arr = np.asarray([v for v in values if v is not None and not np.isnan(v)],
                     dtype=float)

    fig, ax = plt.subplots(figsize=(7.2, 3.9))
    ax.hist(arr, bins=bins, color=SERIES[0], alpha=0.85, lw=0, zorder=2)

    # Stagger the two reference labels vertically - mean and median often sit
    # close enough together that side-by-side text overlaps.
    top = ax.get_ylim()[1]
    for offset, (stat, style, key) in enumerate((
            (np.mean(arr), (0, (4, 3)), "mean"),
            (np.median(arr), (0, (1, 2)), "median"))):
        ax.axvline(stat, color=INK["secondary"], lw=1.2, ls=style, zorder=3)
        ax.annotate(f"{key} {stat:.3f}", xy=(stat, top),
                    xytext=(4, -8 - offset * 12), textcoords="offset points",
                    fontsize=7.5, color=INK["secondary"], va="top")

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Antigens")

    fig.suptitle(title, x=0.012, y=1.02, ha="left", fontsize=12,
                 fontweight="bold", color=INK["primary"])
    fig.text(0.012, 0.95, subtitle or f"n = {arr.size} antigens",
             ha="left", fontsize=8.5, color=INK["muted"])
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fig
