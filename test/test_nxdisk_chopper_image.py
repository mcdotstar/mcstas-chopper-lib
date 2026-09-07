"""What `NXdisk_chopper_image.instr` puts on its imaging plane.

The instrument is a shadowgraph: a flat, collimated, monochromatic field falls on one
parked disc and a PSD directly behind it records what got through. With `abs_out = 1`
the disc body, the hub and the field outside the rim all absorb, so the image *is* the
transmission function -- which makes two things checkable that are otherwise only
checkable by looking at a picture:

  - the lit area is the openings, and nothing else. A component that lost an opening,
    widened one, or drew the annulus between the wrong radii changes the count.
  - a parked disc and a turning one showing the same face give the same image. Those are
    separate branches of the component's TRACE -- `park_angle` against
    `omega * (t - delay)` -- and nothing else in the repository makes them meet.

The source rasters one ray per pixel onto the PSD's own grid, so every assertion here is
on a noiseless image: `N == I` everywhere, and the second test can demand equality pixel
for pixel rather than agreement within some tolerance.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
from mccode_antlr import Flavor

REPO = Path(__file__).resolve().parent.parent
INSTR = REPO / 'NXdisk_chopper_image.instr'

# 256 keeps a run under a second; the geometry does not need more. `-n` must be exactly
# npix * npix or the raster does not cover the frame once.
NPIX = 256
NCOUNT = NPIX * NPIX

# The disc INITIALIZE builds, repeated so the expected area can be worked out from the
# geometry rather than from a previous run. Change one and this file has to change too --
# which is the point: a silent change to the instrument's disc should not pass.
OPENINGS = ((10.0, 18.0), (45.0, 60.0), (100.0, 125.0), (150.0, 190.0), (250.0, 310.0))
RADIUS = 0.35
SLIT_HEIGHT = 0.30
FIELD = 0.72
HUB = RADIUS - SLIT_HEIGHT


def _compiles() -> bool:
    from mccode_antlr.compiler.check import simple_instr_compiles
    return simple_instr_compiles('cc')


requires_compiler = pytest.mark.skipif(
    not _compiles(), reason='these tests run the instrument, so they need a C compiler')


@pytest.fixture(scope='module')
def shadow(tmp_path_factory):
    """The instrument compiled once; call the result for one image.

    Compiling is most of the cost, so the parameter sets share a binary rather than a
    fixture each.
    """
    from mccode_antlr.loader import load_mcstas_instr
    from mccode_antlr.reader.registry import LocalRegistry
    from mccode_antlr.run import mccode_compile, mccode_run_compiled

    # NXdisk_chopper.comp sits beside the instrument; naming the repository explicitly
    # means the tests do not depend on where pytest was started from.
    instr = load_mcstas_instr(INSTR, registries=[LocalRegistry('chopper_lib', str(REPO))])
    work = tmp_path_factory.mktemp('nxdisk_chopper_image')
    binary, target = mccode_compile(instr, work, flavor=Flavor.MCSTAS)
    runs = {}

    def image(**parameters):
        settings = ' '.join(f'{k}={v!r}' if isinstance(v, str) else f'{k}={v}'
                            for k, v in sorted(parameters.items()))
        key = settings
        if key not in runs:
            _, result = mccode_run_compiled(
                binary, target, work / f'run{len(runs)}',
                f'-n {NCOUNT} npix={NPIX} {settings}')
            runs[key] = result['shadow']
        return runs[key]

    return image


def _values(data, field='I'):
    return np.asarray(data[field].values)


def test_the_lit_area_is_the_openings(shadow):
    """The image is annular sectors: five openings, between the hub and the rim."""
    data = shadow(park=0)
    counts = _values(data)

    lit = counts.sum()
    degrees = sum(close - open_ for open_, close in OPENINGS)
    area = degrees / 360.0 * math.pi * (RADIUS ** 2 - HUB ** 2)
    expected = area / FIELD ** 2 * NCOUNT

    # the disc has to actually chop, or "the right area" would say very little
    assert 0.2 * NCOUNT < lit < 0.4 * NCOUNT
    # the remaining difference is the pixels the arcs and the radial edges cut through
    assert lit == pytest.approx(expected, rel=0.01)

    # one unit-weight ray per lit pixel: the raster fell on the monitor's own grid, which
    # is what lets the next test compare images pixel for pixel
    assert _values(data, 'N').sum() == lit
    assert set(np.unique(counts)) == {0.0, 1.0}


def test_a_parked_disc_and_a_turning_one_show_the_same_face(shadow):
    """`park_angle` and `omega * (t - delay)` are separate branches of the same TRACE.

    Rays cross the disc at t = 0, so a turning disc stands at `-360 nu delay` and reaches
    the parked orientation `park` at `delay = -park / (360 nu)`. Both branches then place
    the same openings in front of the same beam, and the images are equal outright -- no
    tolerance, because neither image carries any noise to allow for.
    """
    nu, park = 14.0, 90.0
    parked = _values(shadow(park=park))
    turning = _values(shadow(nu=nu, delay=-park / (360.0 * nu)))

    assert parked.sum() > 0                      # an empty pair would match trivially
    assert np.array_equal(parked, turning)

    # and the equality is a statement about the angle, not about the two branches always
    # agreeing: a disc turned somewhere else does not match
    elsewhere = _values(shadow(nu=nu, delay=-(park + 20.0) / (360.0 * nu)))
    assert not np.array_equal(parked, elsewhere)


# tell pytest to skip the module's tests together when there is nothing to compile with
pytestmark = requires_compiler
