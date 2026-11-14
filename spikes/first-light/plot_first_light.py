#!/usr/bin/env python3
"""Draw the first-light figure (light and dark variants) from the red-noise test output.

    python plot_first_light.py data/pg1302_rednoise.json ../../docs/figures/first-light-pg1302
"""
from __future__ import annotations

import argparse
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator

THEMES = {
    "light": {"surface": "#fcfcfb", "ink": "#0b0b0b", "ink2": "#52514e", "muted": "#898781", "grid": "#e1e0d9",
              "axis": "#c3c2b7", "observed": "#2a78d6", "null": "#eb6834", "errorbar": "#9ec5f4"},
    "dark": {"surface": "#1a1a19", "ink": "#f0efec", "ink2": "#c3c2b7", "muted": "#898781", "grid": "#2c2c2a",
             "axis": "#383835", "observed": "#3987e5", "null": "#d95926", "errorbar": "#184f95"},
}
PX = 0.72  # one CSS pixel in points at 100 px per inch; the figure is saved at twice that


def style(ax, colour):
    ax.set_facecolor(colour["surface"])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(colour["axis"])
        ax.spines[side].set_linewidth(PX)
    ax.tick_params(colors=colour["muted"], labelsize=9, length=3, width=PX)
    ax.grid(True, axis="y", color=colour["grid"], linewidth=PX, linestyle="-")
    ax.set_axisbelow(True)


def render(record, band, out, theme):
    colour = THEMES[theme]
    result, curves = record["bands"][band]["result"], record["bands"][band]["curves"]
    t, mag, err = (np.array(curves[key]) for key in ("t", "mag", "err"))
    period = 1.0 / np.array(curves["freq"])
    power, null_p99 = np.array(curves["power"]), np.array(curves["null_p99"])
    year = 2000.0 + (t - 51544.5) / 365.25
    first, last = (Time(value, format="mjd").strftime("%B %Y") for value in (t.min(), t.max()))

    plt.rcParams["font.family"] = ["DejaVu Sans", "sans-serif"]
    fig, (top, bottom) = plt.subplots(
        2, 1, figsize=(8.0, 6.6), dpi=200, facecolor=colour["surface"],
        gridspec_kw={"hspace": 0.70, "left": 0.085, "right": 0.965, "top": 0.90, "bottom": 0.085})

    style(top, colour)
    top.errorbar(year, mag, yerr=err, fmt="none", ecolor=colour["errorbar"], elinewidth=PX, zorder=2)
    top.plot(year, mag, linestyle="none", marker="o", markersize=8 * PX, markerfacecolor=colour["observed"],
             markeredgecolor=colour["surface"], markeredgewidth=1.5 * PX, zorder=3)
    top.invert_yaxis()
    left, right = np.floor(year.min()), np.ceil(year.max() * 4) / 4 + 0.25
    top.set_xlim(left, right)
    top.xaxis.set_major_locator(FixedLocator([tick for tick in range(int(left), int(right) + 1) if tick <= right]))
    top.set_ylabel(f"{band} magnitude (brighter is up)", color=colour["ink2"], fontsize=9.5)
    top.set_title("PG 1302−102 slowly brightens and fades over seven years of ZTF data",
                  loc="left", color=colour["ink"], fontsize=11.5, fontweight="bold", pad=20)
    top.text(0, 1.045, f"ZTF DR24 {band} band, {result['nights']} nightly means, {first} to {last}",
             transform=top.transAxes, color=colour["ink2"], fontsize=9.5)

    style(bottom, colour)
    bottom.set_xscale("log")
    bottom.plot(period, null_p99, color=colour["null"], linewidth=2 * PX, solid_joinstyle="round",
                solid_capstyle="round", zorder=3,
                label="Level that red noise alone exceeds 1% of the time, at each period")
    bottom.plot(period, power, color=colour["observed"], linewidth=2 * PX, solid_joinstyle="round",
                solid_capstyle="round", zorder=4, label="Observed periodogram")
    claimed = result["claimed_period_days"]
    if claimed:
        bottom.axvline(claimed, color=colour["axis"], linewidth=PX, zorder=1)
        bottom.text(claimed * 0.955, 0.985, f"Published period\n{claimed:,.0f} days", color=colour["ink2"],
                    fontsize=8.5, ha="right", va="top", linespacing=1.35)
    bottom.text(560, 0.42, "Red-noise 1% level", color=colour["ink2"], fontsize=8.5, ha="center", va="bottom")
    bottom.text(1520, 0.30, "Observed", color=colour["ink2"], fontsize=8.5, ha="right", va="center")
    bottom.set_xlim(30, period.max() * 1.02)
    bottom.set_ylim(0, 1.0)
    bottom.xaxis.set_major_locator(FixedLocator([30, 100, 300, 1000, 2500]))
    bottom.xaxis.set_major_formatter(FixedFormatter(["30", "100", "300", "1,000", "2,500"]))
    bottom.xaxis.set_minor_locator(NullLocator())
    bottom.set_xlabel("Trial period (days, log scale)", color=colour["ink2"], fontsize=9.5)
    bottom.set_ylabel("Lomb–Scargle power", color=colour["ink2"], fontsize=9.5)
    bottom.set_title("The strongest peak is not significant once red noise is allowed for",
                     loc="left", color=colour["ink"], fontsize=11.5, fontweight="bold", pad=38)
    exponent = int(round(-np.log10(result["white_noise_fap"])))
    bottom.text(0, 1.045,
                f"A fitted damped random walk (red noise) gives a peak this high in "
                f"{result['red_noise_fap']:.0%} of {result['simulations']:,} simulations.\n"
                f"A standard white-noise test rates the same peak at odds of 1 in $10^{{{exponent}}}$.",
                transform=bottom.transAxes, color=colour["ink2"], fontsize=9.5, va="bottom", linespacing=1.4)
    bottom.legend(loc="upper left", frameon=False, fontsize=9, handlelength=1.6, borderaxespad=0.2,
                  labelcolor=colour["ink2"])

    fig.canvas.draw()
    for text in fig.findobj(matplotlib.text.Text):
        if text.get_text().strip() and text.get_visible():
            box = text.get_window_extent()
            inside = box.x0 >= fig.bbox.x0 and box.x1 <= fig.bbox.x1 and box.y0 >= fig.bbox.y0 and box.y1 <= fig.bbox.y1
            assert inside, f"text falls outside the figure: {text.get_text()!r}"
    fig.savefig(out, facecolor=colour["surface"])
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("results", type=pathlib.Path)
    parser.add_argument("stem", type=pathlib.Path, help="output path without the -light.png / -dark.png suffix")
    parser.add_argument("--band", default="r")
    args = parser.parse_args()
    record = json.loads(args.results.read_text())
    args.stem.parent.mkdir(parents=True, exist_ok=True)
    for theme in THEMES:
        out = args.stem.with_name(f"{args.stem.name}-{theme}.png")
        render(record, args.band, out, theme)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
