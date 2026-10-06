"""What a pencil beam gets past `NXdisk_chopper_transmission.instr`'s one disc.

The image tests check the disc's shape; these check what it does to a beam -- parked or
turning, and where the beam crosses it. A monochromatic pencil beam reaches the disc at
one time, so each run transmits all of its rays or none of them, and most assertions are
exact counts.

Ported from mcstas-readout-master's tests of its `CollectorDiskChopper`, which traced as
this component does and has been retired in its favour.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from mccode_antlr import Flavor

REPO = Path(__file__).resolve().parent.parent
INSTR = REPO / 'NXdisk_chopper_transmission.instr'

RAYS = 200


def _compiles() -> bool:
    from mccode_antlr.compiler.check import simple_instr_compiles
    return simple_instr_compiles('cc')


pytestmark = pytest.mark.skipif(
    not _compiles(), reason='these tests run the instrument, so they need a C compiler')


@pytest.fixture(scope='module')
def transmitted(tmp_path_factory):
    """The instrument compiled once; call the result for the rays one setting passes."""
    from mccode_antlr.loader import load_mcstas_instr
    from mccode_antlr.reader.registry import LocalRegistry
    from mccode_antlr.run import mccode_compile, mccode_run_compiled

    instr = load_mcstas_instr(INSTR, registries=[LocalRegistry('chopper_lib', str(REPO))])
    work = tmp_path_factory.mktemp('nxdisk_chopper_transmission')
    binary, target = mccode_compile(instr, work, flavor=Flavor.MCSTAS)
    runs = []

    def count(**parameters) -> float:
        settings = ' '.join(f'{k}={v}' for k, v in sorted(parameters.items()))
        _, result = mccode_run_compiled(binary, target, work / f'run{len(runs)}',
                                        f'-n {RAYS} {settings}')
        runs.append(settings)
        return float(np.asarray(result['past']['N'].values).sum())

    return count


def test_a_parked_disc_passes_the_beam_through_an_opening(transmitted):
    """The case stock DiskChopper cannot express.

    The disc is stationary with its mark at 180 degrees and a slit spanning it, so the
    opening sits squarely on the beam and rays go through.
    """
    assert transmitted(park=180) == RAYS


def test_a_parked_disc_blocks_when_no_opening_faces_the_beam(transmitted):
    """Same disc, turned so the solid part faces the beam. `DiskChopper` with `nu=0`
    would pass everything here, because it substitutes omega=1e-15 and falls
    permanently open."""
    assert transmitted(park=0) == 0


def test_a_turning_disc_chops_in_time(transmitted):
    """The disc is a clock, not a mask.

    Two runs differing only in `delay`, half a rotation apart, so the same rays meet the
    opening in one and the solid disc in the other. Asserted without saying which is
    which, because that depends on the neutron velocity; what matters is that the two
    are opposite, which a disc that did not chop could not produce.
    """
    counts = [transmitted(nu=14, delay=delay) for delay in (0.0, 1.0 / (2 * 14))]
    assert sorted(counts) == [0, RAYS], counts


def test_the_default_beam_angle_leaves_the_beam_at_the_top(transmitted):
    """Zero is the DiskChopper convention, so nothing already written changes."""
    assert transmitted(edge_open=-10, edge_close=10, beam_angle=0) == RAYS


def test_a_beam_angle_moves_the_beam_round_the_disc(transmitted):
    """A disc hanging above its beam: the opening is at 180, and so is the beam.

    With `beam_angle` the caller no longer turns the whole component about its own z to
    bring that part of the disc to the top -- which is the only thing `DiskChopper` could
    express.
    """
    assert transmitted(beam_angle=180) == RAYS


@pytest.mark.parametrize('label, edges, beam_angle, expected', [
    # the beam moved round to the opening
    ('shifted_beam', (100, 140), 120, RAYS),
    # the opening moved round to the beam: the same disc, said the other way
    ('shifted_edges', (-20, 20), 0, RAYS),
    # and the other direction, which is what fails if the sign is inverted
    ('wrong_way', (100, 140), -120, 0),
])
def test_a_beam_angle_is_a_shift_of_the_openings(transmitted, label, edges, beam_angle,
                                                 expected):
    """What fixes the sense, rather than leaving it to be discovered.

    `beam_angle` enters in the same sense as `slit_edges`, so moving the beam round by B
    is the same disc as moving every edge back by B.
    """
    opening, closing = edges
    assert transmitted(edge_open=opening, edge_close=closing,
                       beam_angle=beam_angle) == expected, label


def test_the_spindle_lies_where_the_angles_put_it(transmitted):
    """Which side of the beam the disc hangs on, tested by what it absorbs.

    An on-axis pencil beam cannot see this -- it is the same distance from a spindle
    above as from one below -- so the disc is moved 40 mm down and `abs_out` is off,
    which makes the two sides behave oppositely: a ray beyond the rim passes, and one
    inside the solid middle does not.

    With the spindle below (`beam_angle=0`) the beam is past the rim and gets through;
    with the disc hanging above it (`beam_angle=180`) the same rays are inside the hub
    and are absorbed. A component that assumed the spindle was always below would pass
    both.

    The beam has width, so it straddles the rim and only some of it clears -- the claim
    is that one side transmits and the other does not, not a precise fraction.
    """
    common = dict(edge_open=0, edge_close=359, abs_out=0, offset_y=-0.04)
    below = transmitted(beam_angle=0, **common)
    above = transmitted(beam_angle=180, **common)
    assert above == 0, (below, above)
    assert below > RAYS / 4, (below, above)


@pytest.mark.parametrize('zero_angle', [0, 90])
def test_zero_angle_does_not_change_what_passes(transmitted, zero_angle):
    """It moves the disc, not the openings.

    Whether a neutron passes is decided on the disc, in the mark's own frame, and
    `zero_angle` says nothing about that -- it only rotates the pickup, and with it the
    spindle, about the component's own axis.
    """
    assert transmitted(park=180, zero_angle=zero_angle) == RAYS
