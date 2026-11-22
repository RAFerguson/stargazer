"""Figures for a blind-search run. Needs matplotlib (`pip install stargazer[figures]`).

Two figures answer the two questions a reader should ask of any null result:

* `calibration_figure`: did the search behave? The share of real objects, and of their decoy twins, at or below
  each probability, against the line a perfectly calibrated null follows.
* `completeness_figure`: could it have found anything? The fraction of planted signals recovered, by shape,
  size and period.

Each is drawn in a light and a dark variant from one palette, with series named on the plot as well as in the
legend so that colour is never the only cue.
"""
from __future__ import annotations

import gzip
import json
import pathlib

import numpy as np

THEMES = {
    "light": {"surface": "#fcfcfb", "ink": "#0b0b0b", "ink2": "#52514e", "muted": "#898781", "grid": "#e1e0d9",
              "axis": "#c3c2b7", "real": "#2a78d6", "decoy": "#eb6834",
              "ramp": ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]},
    "dark": {"surface": "#1a1a19", "ink": "#f0efec", "ink2": "#c3c2b7", "muted": "#898781", "grid": "#2c2c2a",
             "axis": "#383835", "real": "#3987e5", "decoy": "#d95926",
             "ramp": ["#0d366b", "#184f95", "#256abf", "#3987e5", "#6da7ec", "#9ec5f4", "#cde2fb"]},
}
PX = 0.72   # one CSS pixel in points at 100 px per inch; figures are saved at twice that


def load_rows(run_dir: str | pathlib.Path) -> list[dict]:
    run_dir = pathlib.Path(run_dir)
    if (run_dir / "results.json.gz").exists():
        with gzip.open(run_dir / "results.json.gz", "rt") as handle:
            return json.load(handle)
    return json.loads((run_dir / "results.json").read_text())


def _style(ax, colour, grid_axis="y"):
    ax.set_facecolor(colour["surface"])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(colour["axis"])
        ax.spines[side].set_linewidth(PX)
    ax.tick_params(colors=colour["muted"], labelsize=9, length=3, width=PX, which="both")
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=colour["grid"], linewidth=PX, linestyle="-")
    ax.set_axisbelow(True)


def _assert_text_inside(fig) -> None:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    width, height = fig.canvas.get_width_height()
    for text in fig.findobj(lambda artist: hasattr(artist, "get_window_extent") and hasattr(artist, "get_text")):
        if not text.get_text().strip() or not text.get_visible():
            continue
        box = text.get_window_extent(renderer)
        assert box.x0 >= -1 and box.y0 >= -1 and box.x1 <= width + 1 and box.y1 <= height + 1, text.get_text()


def calibration_figure(rows: list[dict], out: str | pathlib.Path, theme: str = "light", title: str = "") -> None:
    """Share of objects at or below each probability, for real objects and decoy twins."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colour = THEMES[theme]
    plt.rcParams["font.family"] = ["DejaVu Sans", "sans-serif"]
    fig, ax = plt.subplots(figsize=(8.0, 4.6), dpi=200, facecolor=colour["surface"])
    fig.subplots_adjust(left=0.095, right=0.80, top=0.80, bottom=0.14)
    _style(ax, colour)
    grid = np.geomspace(1e-3, 1.0, 400)
    ax.plot(grid, grid, color=colour["muted"], linewidth=1.5 * PX, linestyle=(0, (4, 3)), zorder=2)
    ends = {}
    for kind, label in (("decoy", "Decoy twins"), ("real", "Real objects")):
        p = np.sort(np.array([r["p_raw"] for r in rows if r["kind"] == kind and r.get("searched", True)], float))
        if not p.size:
            continue
        share = np.searchsorted(p, grid, side="right") / p.size
        share = np.where(share > 0, share, np.nan)
        ax.plot(grid, share, color=colour[kind], linewidth=2 * PX * 1.4, zorder=3 if kind == "decoy" else 4,
                drawstyle="steps-post", label=f"{label} ({p.size:,})")
        ends[kind] = (label, share[np.isfinite(share)][0], grid[np.isfinite(share)][0])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-3, 1.0)
    ax.set_ylim(5e-4, 1.3)
    ax.set_xlabel("Probability against red noise", color=colour["ink2"], fontsize=10)
    ax.set_ylabel("Share of objects at or below it", color=colour["ink2"], fontsize=10)
    ax.set_xticks([1e-3, 1e-2, 1e-1, 1])
    ax.set_xticklabels(["0.001", "0.01", "0.1", "1"])
    ax.set_yticks([1e-3, 1e-2, 1e-1, 1])
    ax.set_yticklabels(["0.1%", "1%", "10%", "100%"])
    ax.minorticks_off()
    ax.annotate("a calibrated null\nfollows this line", xy=(0.03, 0.03), xytext=(0.06, 0.004), fontsize=9,
                color=colour["ink2"], ha="left", va="center",
                arrowprops={"arrowstyle": "-", "color": colour["muted"], "linewidth": PX})
    legend = ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=9,
                       labelcolor=colour["ink2"], handlelength=1.6, borderaxespad=0.0)
    for handle in legend.legend_handles:
        handle.set_linewidth(2 * PX * 1.4)
    fig.text(0.095, 0.93, title or "Real objects behave like their decoys", color=colour["ink"], fontsize=13,
             fontweight="bold", ha="left", va="center")
    fig.text(0.095, 0.865, "Share of objects whose best period is at least this improbable under red noise",
             color=colour["ink2"], fontsize=9.5, ha="left", va="center")
    _assert_text_inside(fig)
    fig.savefig(out, facecolor=colour["surface"])
    plt.close(fig)


def completeness_figure(cells: list[dict], out: str | pathlib.Path, theme: str = "light", title: str = "") -> None:
    """Fraction of planted signals recovered, one panel per shape, amplitude down and period across."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    colour = THEMES[theme]
    shapes = [s for s in ("sinusoid", "sawtooth", "flare") if any(c["label"] == s for c in cells)]
    amplitudes = sorted({c["amplitude_mag"] for c in cells}, reverse=True)
    periods = sorted({c["period_days"] for c in cells})
    cmap = LinearSegmentedColormap.from_list("ramp", colour["ramp"])
    plt.rcParams["font.family"] = ["DejaVu Sans", "sans-serif"]
    fig, axes = plt.subplots(1, len(shapes), figsize=(8.0, 3.7), dpi=200, facecolor=colour["surface"])
    fig.subplots_adjust(left=0.145, right=0.985, top=0.70, bottom=0.15, wspace=0.10)
    names = {"sinusoid": "Sine wave", "sawtooth": "Sawtooth", "flare": "Narrow pulse"}
    for ax, shape in zip(np.atleast_1d(axes), shapes):
        ax.set_facecolor(colour["surface"])
        for i, amplitude in enumerate(amplitudes):
            for j, period in enumerate(periods):
                cell = next((c for c in cells if c["label"] == shape and c["amplitude_mag"] == amplitude
                             and c["period_days"] == period), None)
                if cell is None or not cell["trials"]:
                    continue
                value = cell["recovered"] / cell["trials"]
                ax.add_patch(plt.Rectangle((j + 0.03, i + 0.03), 0.94, 0.94, facecolor=cmap(value),
                                           edgecolor="none"))
                dark_fill = (value > 0.45) if theme == "light" else (value < 0.55)
                ax.text(j + 0.5, i + 0.42, f"{100 * value:.0f}%", ha="center", va="center", fontsize=11,
                        fontweight="bold", color="#ffffff" if dark_fill else "#0b0b0b")
                ax.text(j + 0.5, i + 0.70, f"{cell['recovered']} of {cell['trials']}", ha="center", va="center",
                        fontsize=7.5, color="#ffffff" if dark_fill else "#0b0b0b")
        ax.set_xlim(0, len(periods))
        ax.set_ylim(len(amplitudes), 0)
        ax.set_xticks(np.arange(len(periods)) + 0.5)
        ax.set_xticklabels([f"{period:.0f} d" for period in periods])
        ax.set_yticks(np.arange(len(amplitudes)) + 0.5)
        ax.set_yticklabels([f"{amplitude:.2f} mag" for amplitude in amplitudes] if ax is np.atleast_1d(axes)[0] else [])
        ax.tick_params(colors=colour["ink2"], labelsize=9, length=0)
        for side in ax.spines.values():
            side.set_visible(False)
        ax.set_title(names[shape], color=colour["ink"], fontsize=10.5, loc="left", pad=6)
        ax.set_xlabel("Period", color=colour["ink2"], fontsize=9)
    np.atleast_1d(axes)[0].set_ylabel("Amplitude", color=colour["ink2"], fontsize=9)
    fig.text(0.03, 0.93, title or "What the search would have found", color=colour["ink"], fontsize=13,
             fontweight="bold", ha="left", va="center")
    fig.text(0.03, 0.855, "Share of signals planted in real light curves that were promoted at the right period",
             color=colour["ink2"], fontsize=9.5, ha="left", va="center")
    _assert_text_inside(fig)
    fig.savefig(out, facecolor=colour["surface"])
    plt.close(fig)


def run_figures(run_dir: str | pathlib.Path, out_stem: str | pathlib.Path, titles: dict | None = None) -> list[str]:
    """Both figures in both themes for one run. Returns the paths written."""
    run_dir = pathlib.Path(run_dir)
    rows = load_rows(run_dir)
    summary = json.loads((run_dir / "summary.json").read_text())
    titles = titles or {}
    written = []
    for theme in THEMES:
        path = f"{out_stem}-calibration-{theme}.png"
        calibration_figure(rows, path, theme, titles.get("calibration", ""))
        written.append(path)
        path = f"{out_stem}-completeness-{theme}.png"
        completeness_figure(summary["completeness"], path, theme, titles.get("completeness", ""))
        written.append(path)
    return written
