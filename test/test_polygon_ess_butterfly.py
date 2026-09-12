"""What `Polygon_ESS_butterfly` emits, and that shaping the source does not change it.

The component refuses to emit outside the region a chopper train transmits, which is only
allowed to make the run *cheaper*, never different. Three ways of running the same
instrument must therefore agree on what reaches the far monitor:

  - `do_region=0`, the source unshaped and every ray emitted;
  - `noise=1`, the region computed but every excluded ray kept anyway;
  - `redraw=1`, excluded rays drawn again from inside the region and every weight
    multiplied by the acceptance.

The third is the one worth testing. It is a different estimator of the same quantity, and
if the acceptance were wrong -- or applied to the wrong rays -- it would agree with neither
of the others.

These run the instrument, so they need a C compiler and skip without one.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest
from mccode_antlr import Flavor

REPO = Path(__file__).resolve().parent.parent
INSTR = REPO / 'Polygon_ESS_butterfly_image.instr'

NCOUNT = 400000


def _compiles() -> bool:
    from mccode_antlr.compiler.check import simple_instr_compiles
    return simple_instr_compiles('cc')


requires_compiler = pytest.mark.skipif(
    not _compiles(), reason='these tests run the instrument, so they need a C compiler')


@pytest.fixture(scope='module')
def butterfly(tmp_path_factory):
    """The instrument compiled once; call the result for one run.

    Compiling is most of the cost, so the parameter sets share a binary.
    """
    from mccode_antlr.loader import load_mcstas_instr
    from mccode_antlr.reader.registry import LocalRegistry
    from mccode_antlr.run import mccode_compile, mccode_run_compiled

    instr = load_mcstas_instr(INSTR, registries=[LocalRegistry('chopper_lib', str(REPO))])
    work = tmp_path_factory.mktemp('polygon_ess_butterfly')
    binary, target = mccode_compile(instr, work, flavor=Flavor.MCSTAS)
    runs = {}

    def run(**parameters):
        settings = ' '.join(f'{k}={v}' for k, v in sorted(parameters.items()))
        if settings not in runs:
            directory = work / f'run{len(runs)}'
            _, result = mccode_run_compiled(
                binary, target, directory, f'-n {NCOUNT} lambda_0=3 {settings}')
            runs[settings] = (result, directory)
        return runs[settings]

    return run


def _monitor(result, name):
    """(intensity, error, counts) from one Monitor_nD of a run.

    McCode writes the three as one `values` header line; reading them from there rather
    than summing the image keeps this the same number the run printed.
    """
    intensity, error, counts = result[name].metadata['values'].split()
    return float(intensity), float(error), float(counts)


def _region(directory):
    """The region the component wrote beside its data."""
    with open(Path(directory) / 'source.json') as file:
        return json.load(file)


@requires_compiler
def test_the_region_is_written_and_self_consistent(butterfly):
    _, directory = butterfly(redraw=0)
    region = _region(directory)

    assert region['polygons'], 'the train transmits something, so there is a polygon'
    assert region['transmitted_area'] > 0
    assert region['acceptance'] == pytest.approx(
        region['transmitted_area'] / region['sampled']['area'], rel=1e-12)
    assert sum(p['area'] for p in region['polygons']) == pytest.approx(
        region['transmitted_area'], rel=1e-12)
    for polygon in region['polygons']:
        assert len(polygon['vertices']) >= 3
        for inverse_velocity, time in polygon['vertices']:
            assert region['sampled']['inverse_velocity'][0] - 1e-12 <= inverse_velocity
            assert inverse_velocity <= region['sampled']['inverse_velocity'][1] + 1e-12
            assert region['sampled']['time'][0] - 1e-12 <= time
            assert time <= region['sampled']['time'][1] + 1e-12


@requires_compiler
def test_the_acceptance_is_what_the_run_measures(butterfly):
    """The geometry computes it; the run counts it. The two are independent."""
    result, directory = butterfly(redraw=0)
    region = _region(directory)
    _, _, emitted = _monitor(result, 'emission')
    # Every ray the source drew is either emitted or absorbed for being outside the
    # region, so the emitted fraction of ncount is the acceptance measured from the run.
    measured = emitted / NCOUNT
    assert measured == pytest.approx(region['acceptance'], rel=0.02)


@requires_compiler
@pytest.mark.parametrize('shaped', ['noise=1', 'redraw=1'])
def test_shaping_the_source_does_not_change_what_gets_through(butterfly, shaped):
    """Whatever the source does with the region, the discs pass the same beam."""
    key, value = shaped.split('=')
    plain, _ = butterfly(do_region=0)
    other, _ = butterfly(**{key: int(value)})

    a, a_error, _ = _monitor(plain, 'transmitted')
    b, b_error, _ = _monitor(other, 'transmitted')
    combined = math.hypot(a_error, b_error)
    assert abs(a - b) < 4 * combined, f'{a:.6g} +- {a_error:.3g} vs {b:.6g} +- {b_error:.3g}'


@requires_compiler
def test_redrawing_buys_statistics(butterfly):
    """The whole point: the same answer out of the same ncount, with more counts in it.

    Resampling spends no ray on a region the discs would have shut, so the surviving count
    rises by roughly one over the acceptance.
    """
    plain, _ = butterfly(do_region=0)
    redrawn, directory = butterfly(redraw=1)
    acceptance = _region(directory)['acceptance']

    _, plain_error, plain_counts = _monitor(plain, 'transmitted')
    _, redrawn_error, redrawn_counts = _monitor(redrawn, 'transmitted')

    assert redrawn_counts == pytest.approx(plain_counts / acceptance, rel=0.1)
    # more counts is less error, by the square root of the ratio
    assert redrawn_error < plain_error * math.sqrt(1.5 * acceptance)


@requires_compiler
def test_a_path_spread_only_widens_the_region(butterfly):
    """Extra flight path can only let more through, and only in time.

    The bands it opens are wider in inverse velocity too -- a slower neutron is smeared
    more -- but nothing is ever removed, so the area cannot fall.
    """
    _, straight = butterfly(redraw=0)
    _, spread = butterfly(spread=1e-3)
    tight = _region(straight)
    loose = _region(spread)
    assert loose['transmitted_area'] >= tight['transmitted_area']
    assert loose['acceptance'] >= tight['acceptance']

    low, high = tight['inverse_velocity_bands'][0]
    spread_low, spread_high = loose['inverse_velocity_bands'][0]
    assert spread_low <= low + 1e-12
    assert spread_high >= high - 1e-12
