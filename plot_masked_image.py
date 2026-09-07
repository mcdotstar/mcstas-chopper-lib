#!/usr/bin/env python
"""Four panels of one Masked_ESS_butterfly_image run, on one set of axes.

    mcstas-antlr Masked_ESS_butterfly_image.instr
    ./Masked_ESS_butterfly_image.out -n 2000000 lambda_0=3 -d run
    python plot_masked_image.py run -o mask.png

The run leaves four pictures of the same (wavelength, emission time) plane, and
they only mean anything against each other:

    source.total      every sample the source drew, on the mask's own grid
    source.mask       which of those bins chopper-lib said a neutron could pass
    emission.L_U1     what the source emitted, on the monitor's finer grid
    transmitted.L_U1  what the discs then passed, on those same monitor axes

Read left to right, top to bottom, that is the whole of what the component does:
the source samples a frame, the library cuts a shape out of it, the source emits
only inside that shape, and the discs the shape was computed from pass most --
not all -- of what came out. The last panel differs from the third because a
disc's opening is angular and the mask's is not; see the instrument's own notes.

The monitor panels plot intensity -- the summed ray weight, which is what the
source's own `.total` file holds and the only channel comparable with it. `--counts`
plots the ray count instead, which answers a different question: how well sampled a
bin is, rather than how much of the beam is in it. The two look alike here only
because the source samples wavelength and emission time uniformly.

Panels one, three and four are magnitudes and share one light-to-dark blue ramp,
so darker is always more; empty bins are left the colour of the paper rather than
the pale end of the ramp. The emission and transmission panels share a scale as
well as a ramp, which is what makes the discs' bite out of the ribbon legible.
The mask is a state rather than a magnitude, so it takes the two ends of the same
ramp with its colourbar labelled in words.

Requires mccode-antlr with scipp for the monitor files; the mask and total files
are plain text this reads directly.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

# Steps of one blue ramp, light to dark. Sequential data gets one hue and a
# monotone lightness; a rainbow would put false edges in a smooth field.
RAMP = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b']

LIGHT = {'surface': '#fcfcfb', 'ink': '#0b0b0b', 'muted': '#52514e', 'line': '#c9c8c4'}
DARK = {'surface': '#1a1a19', 'ink': '#ffffff', 'muted': '#c3c2b7', 'line': '#4a4a47'}

V2K = 1.58825361e-3     # v [m/s] to k [1/AA], as the McCode runtime defines it


def wavelength_of(inverse_velocity):
    """Wavelength in AA of an inverse velocity in s/m."""
    return np.asarray(inverse_velocity) / (V2K / 2 / np.pi)


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------

def read_chopper_grid(path: Path):
    """A chopper-lib `.mask` or `.total` file.

    Both are written by the same pair of functions and share a header: two axes
    as bin edges, then one row per time bin. The values are the only difference
    -- a mask holds 0 or 1, a total holds the summed ray weight.
    """
    lines = [line.rstrip('\n') for line in path.read_text().splitlines()]
    edges = [np.array([float(x) for x in line.lstrip('#').split()])
             for line in lines if line.startswith('#') and line.lstrip('#')[:1].isdigit()]
    if len(edges) != 2:
        raise SystemExit(f'{path} does not carry the two axes a chopper-lib grid file has')
    inverse_velocity, time = edges
    values = np.array([[float(x) for x in line.split()]
                       for line in lines if line and not line.startswith('#')])
    if values.shape != (time.size - 1, inverse_velocity.size - 1):
        raise SystemExit(f'{path}: {values.shape} values do not fit its axes')
    # rows are time bins and columns inverse velocity bins, which is already
    # (y, x) for an image of wavelength against time
    return values, wavelength_of(inverse_velocity), time


def read_monitor(path: Path, channel: str = 'I'):
    """A Monitor_nD file, as (values, wavelength edges, time edges, parameters).

    *channel* is `'I'` for the intensity -- the summed ray weight, which is what the
    source's own `.total` file holds and the only channel comparable with it -- or
    `'N'` for the ray count, which says how well sampled a bin is rather than how
    much of the beam is in it.
    """
    from mccode_antlr.loader.datfile import read_mccode_dat

    dat = read_mccode_dat(path)
    counts = np.asarray(dat.dataset[channel].values, dtype=float)
    lo_l, hi_l, lo_t, hi_t = [float(x) for x in dat.metadata['xylimits'].split()]
    n_t, n_l = counts.shape
    return (counts,
            np.linspace(lo_l, hi_l, n_l + 1),
            np.linspace(lo_t, hi_t, n_t + 1),
            dat.parameters)


def find(directory: Path, pattern: str, what: str) -> Path:
    """The one file in *directory* matching *pattern*."""
    found = sorted(directory.glob(pattern))
    if not found:
        raise SystemExit(f'No {what} ({pattern}) in {directory}')
    return found[0]


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------

def ramp(theme):
    """The blue ramp as a colour map, with empty bins painted as the surface.

    On a dark surface the ramp is not simply flipped: it runs from the surface
    itself up through the same blues to their palest step, so that a low value
    reads as empty ground in either theme rather than as a field of deep blue.
    """
    from matplotlib.colors import LinearSegmentedColormap

    steps = RAMP if theme is LIGHT else [theme['surface'], *RAMP[::-1]]
    cmap = LinearSegmentedColormap.from_list('chopper', steps)
    cmap = cmap.copy()
    cmap.set_under(theme['surface'])
    return cmap


def two_state(theme):
    """The two ends of the ramp, for a field that is only ever included or not."""
    from matplotlib.colors import ListedColormap

    return ListedColormap([theme['surface'], RAMP[-1] if theme is LIGHT else RAMP[0]])


def draw(ax, values, lam_edges, t_edges, cmap, theme, vmin, vmax, title, subtitle):
    """One panel: an image in (wavelength, emission time), and what it is."""
    image = ax.pcolormesh(
        lam_edges, t_edges * 1e3, values,
        cmap=cmap, vmin=vmin, vmax=vmax, shading='flat', rasterized=True,
    )
    ax.set_title(title, color=theme['ink'], fontsize=10.5, loc='left', pad=6)
    # the file this panel came from, inside the frame where it cannot collide with
    # the panel above it, and boxed so it stays legible over a busy corner
    ax.text(0.975, 0.955, subtitle, transform=ax.transAxes, color=theme['muted'],
            fontsize=8, family='monospace', ha='right', va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor=theme['surface'],
                      edgecolor=theme['line'], linewidth=0.6, alpha=0.9))
    ax.set_facecolor(theme['surface'])
    ax.tick_params(colors=theme['muted'], labelsize=8.5, length=3)
    for spine in ax.spines.values():
        spine.set_color(theme['line'])
        spine.set_linewidth(0.8)
    return image


def bar(fig, image, axes, theme, label, ticks=None, ticklabels=None, shrink=1.0):
    """A colourbar for one panel or a pair of them, labelled with its units."""
    cb = fig.colorbar(image, ax=axes, fraction=0.055, pad=0.02, ticks=ticks,
                      shrink=shrink)
    cb.set_label(label, color=theme['muted'], fontsize=8.5)
    cb.ax.tick_params(colors=theme['muted'], labelsize=8, length=3)
    cb.outline.set_edgecolor(theme['line'])
    cb.outline.set_linewidth(0.8)
    if ticklabels is not None:
        cb.ax.set_yticklabels(ticklabels)
    return cb


def caption(parameters: dict) -> str:
    """The train, read back from the parameters the monitor file carries."""
    def value(name, default=None):
        raw = parameters.get(name, default)
        try:
            return float(raw)
        except (TypeError, ValueError):
            return default

    near, far = value('near_path'), value('far_path')
    if near is None or far is None:
        return ''
    return (f"discs at {near:g} m ({value('near_width'):g}°) and {far:g} m "
            f"({value('far_width'):g}°) turning at {value('near_nu'):g} Hz, "
            f"set for {value('lambda_0'):g} Å leaving at "
            f"{1e3 * value('t_0'):g} ms")


def interesting_box(panels, margin=0.25):
    """The corner of the frame the run actually used, with room around it.

    Every panel but the sampling one is empty over most of the plane -- the frame
    is the pulse and the whole wavelength band, and the ribbon is a small part of
    it. Plotting the whole frame would leave four mostly blank squares.
    """
    lam_lo, lam_hi, t_lo, t_hi = np.inf, -np.inf, np.inf, -np.inf
    for values, lam_edges, t_edges in panels:
        rows, cols = np.nonzero(values > 0)
        if not rows.size:
            continue
        lam_lo = min(lam_lo, lam_edges[cols.min()])
        lam_hi = max(lam_hi, lam_edges[cols.max() + 1])
        t_lo = min(t_lo, t_edges[rows.min()])
        t_hi = max(t_hi, t_edges[rows.max() + 1])
    if not np.isfinite([lam_lo, lam_hi, t_lo, t_hi]).all():
        return None
    pad_l, pad_t = margin * (lam_hi - lam_lo), margin * (t_hi - t_lo)
    return (lam_lo - pad_l, lam_hi + pad_l, 1e3 * (t_lo - pad_t), 1e3 * (t_hi + pad_t))


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__.split('\n\n')[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('directory', type=Path, nargs='?', default=Path('.'),
                        help='the run directory holding all four files')
    parser.add_argument('-o', '--output', type=Path, default=Path('masked_image.png'),
                        help='image to write; any format matplotlib knows')
    parser.add_argument('--name', default='source',
                        help='the source component name, which names its two files')
    parser.add_argument('--emission', default='emission*',
                        help='the near monitor file, or a pattern for it')
    parser.add_argument('--transmitted', default='transmitted*',
                        help='the far monitor file, or a pattern for it')
    parser.add_argument('--dpi', type=int, default=140, help='output resolution')
    parser.add_argument('--dark', action='store_true',
                        help='draw on a dark surface, with the ramp reversed to suit')
    parser.add_argument('--counts', action='store_true',
                        help='plot the monitors\' ray counts rather than their intensity')
    parser.add_argument('--zoom', action='store_true',
                        help='crop to the part of the frame the run used, with a margin')
    parser.add_argument('--title', default=None, help='override the figure title')
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    root = args.directory
    total_file = find(root, f'{args.name}.total', 'total file')
    mask_file = find(root, f'{args.name}.mask', 'mask file')
    emission_file = find(root, args.emission, 'near monitor file')
    transmitted_file = find(root, args.transmitted, 'far monitor file')

    total, total_lam, total_t = read_chopper_grid(total_file)
    mask, mask_lam, mask_t = read_chopper_grid(mask_file)
    channel = 'N' if args.counts else 'I'
    emitted, em_lam, em_t, parameters = read_monitor(emission_file, channel)
    passed, tr_lam, tr_t, _ = read_monitor(transmitted_file, channel)

    theme = DARK if args.dark else LIGHT
    cmap = ramp(theme)

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.6), dpi=args.dpi,
                             sharex=True, sharey=True, facecolor=theme['surface'],
                             layout='constrained')
    (ax_total, ax_mask), (ax_emitted, ax_passed) = axes

    # A shared scale for the two monitor panels: they count the same thing on the
    # same axes, so the discs' bite out of the ribbon has to be a change in colour
    # rather than a change in what the colours mean.
    counts_max = max(emitted.max(), passed.max())
    # `vmin` a hair above zero, so empty bins take the map's `under` colour and
    # read as paper rather than as the pale end of the ramp
    eps = 1e-9

    image_total = draw(ax_total, total, total_lam, total_t, cmap, theme,
                       eps, total.max(), 'Every sample the source drew',
                       total_file.name)
    image_mask = draw(ax_mask, mask.astype(float), mask_lam, mask_t, two_state(theme),
                      theme, -0.5, 1.5, 'What chopper-lib said could pass',
                      mask_file.name)
    image_emitted = draw(ax_emitted, emitted, em_lam, em_t, cmap, theme,
                         eps, counts_max, 'What the source emitted',
                         emission_file.name)
    draw(ax_passed, passed, tr_lam, tr_t, cmap, theme,
         eps, counts_max, 'What the discs then passed',
         transmitted_file.name)

    bar(fig, image_total, ax_total, theme, 'summed ray weight')
    # two states, so a short bar that reads as a key rather than a scale
    bar(fig, image_mask, ax_mask, theme, '', ticks=[0, 1],
        ticklabels=['excluded', 'included'], shrink=0.45)
    # one bar for the pair: they count the same thing, so the discs' bite has to
    # show as a change of colour and not a change of scale
    bar(fig, image_emitted, [ax_emitted, ax_passed], theme,
        ('rays per bin, one scale' if args.counts else 'intensity per bin, one scale'))

    for ax in axes[1]:
        ax.set_xlabel('Wavelength [Å]', color=theme['muted'], fontsize=9.5)
    for ax in axes[:, 0]:
        ax.set_ylabel('Emission time [ms]', color=theme['muted'], fontsize=9.5)

    if args.zoom:
        box = interesting_box([(mask, mask_lam, mask_t),
                               (emitted, em_lam, em_t),
                               (passed, tr_lam, tr_t)])
        if box is not None:
            ax_total.set_xlim(box[0], box[1])
            ax_total.set_ylim(box[2], box[3])

    # the layout engine packs the axes into this rectangle, leaving the strip above
    # it for the title and the line under it
    fig.get_layout_engine().set(rect=(0, 0, 1, 0.915))
    fig.text(0.012, 0.985, args.title or 'Masked_ESS_butterfly_image',
             color=theme['ink'], fontsize=13, ha='left', va='top')
    line = caption(parameters)
    if line:
        fig.text(0.012, 0.945, line, color=theme['muted'], fontsize=9,
                 ha='left', va='top')

    fig.savefig(args.output, facecolor=theme['surface'])
    plt.close(fig)
    print(f'Wrote {args.output}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
