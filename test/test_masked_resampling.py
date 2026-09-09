"""What `Masked_ESS_butterfly(resample=1)` measures, and what it refuses to measure.

A chopper mask tells the source which (inverse velocity, emission time) cells the discs
can pass. Absorbing everything outside them is right but expensive: a train passing a few
percent of the plane spends the rest of `ncount` on rays that die where they are born.
`resample=1` draws those rays again from inside the allowed region instead, and pays for
the restriction by multiplying every ray weight by the region's share of the plane.

That factor is the whole claim, and a wrong one is invisible in a picture -- the image
looks right and every number in it is off by a constant. So the tests here are about
intensity, not shape:

  - the three ways of running the same instrument (no mask, mask-and-absorb, mask-and-
    redraw) measure the same transmitted intensity, while the third gets there with a
    much smaller error bar. Agreement alone would pass with the mask switched off, and a
    smaller error bar alone would pass with the weight factor missing; the pair will not.
  - nothing at all is emitted outside the mask when redrawing, which is what makes the
    factor a restriction rather than a scaling.
  - the acceptance the component computes from the mask's geometry matches the fraction
    of draws that actually landed in it, which the component measures separately and
    reports. Neither number is derived from the other.

and then about the two configurations where the argument does not hold -- time focusing
and more than one pulse -- which INITIALIZE has to refuse rather than silently bias.
"""
from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
from mccode_antlr import Flavor

REPO = Path(__file__).resolve().parent.parent
IMAGE = REPO / 'Masked_ESS_butterfly_image.instr'

# Enough rays for the error bars to mean something and still under a few seconds. It also
# has to divide by every rank count the MPI tests below use: mccode_main hands each node
# floor(ncount / nodes) rays, so an indivisible count drops the remainder and the draw total
# stops being exactly NCOUNT.
NCOUNT = 400_000
MPI_RANKS = (1, 2, 4)
# Fixed, so a failure is reproducible and a marginal run does not flap in CI.
SEED = 20250908
# The wavelength the default two-disc train is set for; anything the mask passes will do.
LAMBDA_0 = 3.0


def _compiles() -> bool:
    from mccode_antlr.compiler.check import simple_instr_compiles
    return simple_instr_compiles('cc')


requires_compiler = pytest.mark.skipif(
    not _compiles(), reason='these tests run the instrument, so they need a C compiler')


def _compile(instr_path: Path, work: Path, ranks: int | None = None):
    """The instrument built once, and a way to run it that keeps its stdout.

    `Masked_ESS_butterfly` reports its acceptance on stdout and refuses on stderr, and
    both are part of what is under test here, so the runs go through the binary directly
    rather than through `mccode_run_compiled`.
    """
    from mccode_antlr.loader import load_mcstas_instr
    from mccode_antlr.reader.registry import LocalRegistry
    from mccode_antlr.run import mccode_compile

    instr = load_mcstas_instr(instr_path, registries=[LocalRegistry('chopper_lib', str(REPO))])
    target = None if ranks is None else {'mpi': True, 'count': ranks}
    binary, _ = mccode_compile(instr, work, flavor=Flavor.MCSTAS, target=target)
    runs = {}

    def run(ranks=None, **parameters):
        settings = [f'{k}={v}' for k, v in sorted(parameters.items())]
        key = f'{ranks} ' + ' '.join(settings)
        if key not in runs:
            directory = work / f'run{len(runs)}'
            # mpirun is invoked here rather than through mccode_run_compiled so a refused run
            # comes back as a return code to assert on, the way the serial runs above do,
            # instead of as a RuntimeError.
            launch = [] if ranks is None else ['mpirun', '-np', str(ranks)]
            runs[key] = subprocess.run(
                [*launch, str(binary), '-n', str(NCOUNT), '-s', str(SEED),
                 '-d', str(directory), *settings],
                capture_output=True, text=True)
            runs[key].directory = directory
        return runs[key]

    return run


@pytest.fixture(scope='module')
def image(tmp_path_factory):
    return _compile(IMAGE, tmp_path_factory.mktemp('masked_resampling'))


def _monitor(completed, name):
    """The total intensity, its error and the count in one of the instrument's images."""
    from mccode_antlr.loader import read_mccode_dat
    assert completed.returncode == 0, completed.stderr
    data = read_mccode_dat(completed.directory / f'{name}.L_U1')
    intensity = float(data['I'].values.sum())
    # Bin errors are independent, so the total's error is their sum in quadrature; the
    # dat file carries them as variances, which is that sum already.
    error = float(np.sqrt(data['I'].variances.sum()))
    return intensity, error, float(data['N'].values.sum())


def _psd(completed):
    """The total intensity and count on the guard instrument's monitor."""
    from mccode_antlr.loader import read_mccode_dat
    assert completed.returncode == 0, completed.stderr
    data = read_mccode_dat(completed.directory / 'psd.dat')
    return float(data['I'].values.sum()), float(data['N'].values.sum())


def _grid(completed, extension):
    """One of the component's own (inverse velocity, time) files as an array."""
    return np.loadtxt(completed.directory / f'source{extension}', comments='#')


@requires_compiler
def test_redrawing_measures_what_absorbing_measures(image):
    """The same transmitted intensity three ways, and a smaller error bar for the third.

    `do_mask=0` emits across the whole sampled plane and lets the discs do the cutting,
    which is the answer with no mask involved at all. `do_mask=1` is the mask applied by
    absorbing. `redraw=1` is the mask applied by resampling. All three are estimates of
    one number, so they have to agree to within their own error bars -- and the third's
    error bar has to be the small one, or the resampling bought nothing.
    """
    unmasked = _monitor(image(lambda_0=LAMBDA_0, do_mask=0), 'transmitted')
    absorbed = _monitor(image(lambda_0=LAMBDA_0, do_mask=1), 'transmitted')
    redrawn = _monitor(image(lambda_0=LAMBDA_0, redraw=1), 'transmitted')

    # The discs must actually be cutting something, or "they agree" says very little
    assert 0 < redrawn[0]
    assert unmasked[2] < NCOUNT / 2

    for other in (unmasked, absorbed):
        combined = np.hypot(other[1], redrawn[1])
        assert abs(redrawn[0] - other[0]) < 4 * combined, (redrawn, other)

    # Where the resampling actually pays: the same ncount lands far more rays past the
    # discs, and the error bar falls with them.
    assert redrawn[2] > 2 * absorbed[2]
    assert redrawn[1] < 0.75 * absorbed[1]


@requires_compiler
def test_redrawing_emits_only_inside_the_mask(image):
    """`resample=1` is a restriction of the sampling, not a rescaling of it.

    The component writes what it emitted on the mask's own grid, so this is checkable
    cell by cell. The `do_mask=0` half is not decoration: without it an `.emitted` file
    that happened to be zero everywhere outside the mask for some other reason would pass.
    """
    mask = _grid(image(lambda_0=LAMBDA_0, redraw=1), '.mask')
    emitted = _grid(image(lambda_0=LAMBDA_0, redraw=1), '.emitted')
    assert mask.shape == emitted.shape

    assert np.all(emitted[mask == 0] == 0.0)
    assert np.count_nonzero(emitted[mask == 1]) > 0

    # and the mask is a real constraint on this instrument, not the whole frame
    assert 0 < np.count_nonzero(mask == 0)
    spilled = _grid(image(lambda_0=LAMBDA_0, do_mask=0), '.emitted')
    assert np.count_nonzero(spilled[mask == 0]) > 0


@requires_compiler
def test_the_acceptance_is_the_fraction_of_draws_it_claims(image):
    """The weight factor against the sampling it is supposed to describe.

    The acceptance is computed in INITIALIZE from the mask's geometry; the sampled figure
    is counted during the run, unweighted, from where the source's *first* draw landed
    before any redrawing. Nothing connects the two but the claim under test, so they
    agree only if that claim is right.
    """
    completed = image(lambda_0=LAMBDA_0, redraw=1)
    match = re.search(r'mask acceptance ([\d.eE+-]+); sampled ([\d.eE+-]+)', completed.stdout)
    assert match is not None, completed.stdout
    acceptance, sampled = float(match.group(1)), float(match.group(2))

    # A binomial fraction from NCOUNT draws; three sigma either way.
    sigma = np.sqrt(acceptance * (1 - acceptance) / NCOUNT)
    assert abs(sampled - acceptance) < 3 * sigma, (acceptance, sampled, sigma)
    # and the mask really is a restriction, so this is not two ones being compared
    assert 0.0 < acceptance < 1.0


# The window sweeps a band `tfocus_dist * inverse_velocity_range` longer than itself and
# reaches back before the pulse starts, so this drives both the widening and the clip.
TIME_FOCUS = dict(tf_width=1e-3, tf_time=0.004, tf_dist=5.0)

# A minimal instrument, because the imaging one exposes neither time focusing nor the
# pulse count and has no reason to. One wide disc, so the mask allows plenty and nothing
# but the guard under test can stop the run.
GUARD_INSTR = """
DEFINE INSTRUMENT resample_guard(tf_width=0, tf_time=0, tf_dist=0, int pulses=1,
                                 int redraw=1, noise=0, int use=1, tmax=3)
DECLARE %{
chopper_parameters * train;
double * train_as_doubles;
double edges[2];
%}
INITIALIZE %{
edges[0] = -30.0; edges[1] = 30.0;
train = (chopper_parameters *) calloc(1, sizeof(chopper_parameters));
train[0] = (chopper_parameters){14.0, 0.004, 0.0, 2, edges, 5.0, 0.0};
train_as_doubles = (double *) train;
%}
TRACE
COMPONENT origin = Progress_bar() AT (0, 0, 0) ABSOLUTE
COMPONENT source = Masked_ESS_butterfly(
  sector="N", beamline=1, Lmin=1.0, Lmax=5.0, yheight=0.03, cold_frac=0.5,
  dist=5.0, focus_xw=0.1, focus_yh=0.1, n_pulses=pulses,
  tfocus_dist=tf_dist, tfocus_time=tf_time, tfocus_width=tf_width,
  choppers=train_as_doubles, chopper_count=1,
  inverse_velocity_bin=1e-5, time_bin=1e-4,
  tmax_multiplier=tmax,
  filename="guard", noise_fraction=noise, use_mask=use, resample=redraw
) AT (0, 0, 0) ABSOLUTE
COMPONENT psd = PSD_monitor(xwidth=0.5, yheight=0.5, nx=10, ny=10, filename="psd",
                            restore_neutron=1)
AT (0, 0, 1) RELATIVE source
FINALLY
%{
if (train) free(train);
%}
END
"""


@pytest.fixture(scope='module')
def guard(tmp_path_factory):
    work = tmp_path_factory.mktemp('resample_guard')
    path = work / 'resample_guard.instr'
    path.write_text(GUARD_INSTR)
    return _compile(path, work)


@requires_compiler
def test_the_guarded_instrument_runs_when_nothing_is_wrong(guard):
    """The baseline the refusals are refusals *from*."""
    completed = guard(redraw=1)
    assert completed.returncode == 0, completed.stderr
    assert 'acceptance' in completed.stdout


@requires_compiler
@pytest.mark.parametrize('parameters,expected', [
    # The emission time is drawn in a window centred on tfocus_time - tfocus_dist/vz, so
    # it is not independent of the velocity and no constant factor puts it back.
    (dict(TIME_FOCUS), 'time focusing'),
    # Resampling emits nothing outside the mask, so there is no leak to set a rate for.
    (dict(noise=0.5), 'contradictory'),
])
def test_resampling_refuses_what_it_cannot_correct(guard, parameters, expected):
    completed = guard(redraw=1, **parameters)
    assert completed.returncode != 0, completed.stdout
    assert expected in completed.stderr, completed.stderr


@requires_compiler
@pytest.mark.parametrize('parameters', [dict(TIME_FOCUS), dict(noise=0.5)])
def test_the_refusals_are_about_resampling_and_nothing_else(guard, parameters):
    """Each configuration refused above is fine with the mask applied by absorbing.

    Without this the parametrised test above would also pass if the component had simply
    stopped working for those parameters.
    """
    completed = guard(redraw=0, **parameters)
    assert completed.returncode == 0, completed.stderr


@requires_compiler
def test_time_focusing_leaves_a_beam_to_mask(guard):
    """The grid has to span the band the focusing window sweeps, not the window.

    The window is `tfocus_width` wide but its centre slides with `tfocus_dist / vz`, so the
    band it covers is `tfocus_dist * inverse_velocity_range` longer -- and it reaches back
    before the pulse begins, so it also has to be clipped to what the source can emit.
    Sizing the grid by `tfocus_width` alone puts it around the slowest neutron's window and
    nowhere near anything else, and TRACE's bounds check then absorbs every ray whether the
    mask is switched on or not.

    So what is checked here is not that the mask passes something; it is that the *unmasked*
    run passes something, which is a statement about the grid and nothing else.
    """
    unmasked = _psd(guard(redraw=0, use=0, **TIME_FOCUS))
    masked = _psd(guard(redraw=0, use=1, **TIME_FOCUS))

    assert unmasked[1] > 0, 'the grid absorbed the whole beam before the mask saw it'
    # a decent share of the pulse survives: focusing rejects the neutrons whose window falls
    # before t = 0, and nothing else should be lost
    assert unmasked[1] > NCOUNT / 4
    # the mask may cut into that, but it cannot be the thing emptying it
    assert masked[1] > 0
    assert masked[0] <= unmasked[0] * (1 + 1e-9)


@requires_compiler
@pytest.mark.parametrize('redraw', [0, 1])
def test_more_pulses_do_not_lose_intensity(guard, redraw):
    """The pulse offset is taken back off before the mask is consulted.

    ESS_butterfly picks a pulse per ray and adds its offset to `t`; a mask indexed on that
    sees a time outside its grid for every ray but the first pulse's, and TRACE's bounds
    check absorbs them -- so the intensity comes out low by a factor of `n_pulses`.

    Spreading rays over pulses does not change what the source emits in total, and
    `floor(n_pulses * rand01())` is drawn whether there is one pulse or several, so the two
    runs see identical random streams and identical weights. The intensities are therefore
    equal outright rather than equal within statistics, which is a far sharper thing to
    assert than a tolerance would be.
    """
    one = _psd(guard(redraw=redraw, pulses=1))
    three = _psd(guard(redraw=redraw, pulses=3))

    assert one[1] > 0
    assert three[0] == pytest.approx(one[0], rel=1e-12)
    assert three[1] == one[1]


@requires_compiler
def test_a_pulse_offset_that_cannot_be_recovered_is_refused(guard):
    """Splitting `t` needs the emission window to be shorter than the gap between pulses.

    It is, by a factor of eight, at the default `tmax_multiplier` of 3 -- but the parameter
    is free, and at a large enough value the pulses overlap and no arithmetic separates the
    offset from the time it was added to. Guessing there would put rays in the wrong bins
    silently, so INITIALIZE stops instead.
    """
    # 1 / 14 Hz is 71.4 ms; 30 * 2.857 ms is 85.7 ms, which no longer fits
    completed = guard(redraw=0, pulses=3, tmax=30)
    assert completed.returncode != 0, completed.stdout
    assert 'pulse period' in completed.stderr, completed.stderr

    # and it is the overlap that is refused, not the multiplier: one pulse is fine
    assert guard(redraw=0, pulses=1, tmax=30).returncode == 0


# ---------------------------------------------------------------------------
# MPI
#
# Every node runs SAVE -- mccode_main calls finally() on all of them and finally() calls
# save() unconditionally -- and mcuse_dir hands them all one output directory. So the
# component has to sum its three grids across the nodes and let master alone write them, and
# the two halves fail in different ways: without the gate each file holds one interleaved
# copy per node, and without the sum each copy holds one node's share.
# ---------------------------------------------------------------------------
def _mpi_available() -> bool:
    from mccode_antlr.compiler.check import simple_instr_compiles
    return shutil.which('mpirun') is not None and simple_instr_compiles('mpi/cc')


requires_mpi = pytest.mark.skipif(
    not _mpi_available(), reason='these tests need an MPI compiler and mpirun')


@pytest.fixture(scope='module')
def parallel(tmp_path_factory):
    """The guard instrument built against mpicc, run under `mpirun -np ranks`."""
    work = tmp_path_factory.mktemp('masked_mpi')
    path = work / 'resample_guard.instr'
    path.write_text(GUARD_INSTR)
    return _compile(path, work, ranks=max(MPI_RANKS))


def _draws(completed):
    """The whole run's draw count, as the component reports it after reducing."""
    assert completed.returncode == 0, completed.stderr
    match = re.search(r'over ([\d.eE+-]+) draws', completed.stdout)
    assert match is not None, completed.stdout
    return float(match.group(1))


def _headers(completed, extension):
    """How many times the component wrote one of its files into the run directory."""
    assert completed.returncode == 0, completed.stderr
    text = (completed.directory / f'guard{extension}').read_text()
    return text.count('# Chopper mask file generated by chopper-lib')


@requires_mpi
@pytest.mark.parametrize('ranks', MPI_RANKS)
@pytest.mark.parametrize('extension', ['.mask', '.total', '.emitted'])
def test_each_output_file_is_written_once(parallel, ranks, extension):
    """Master writes; the other nodes do not.

    Every node reaches SAVE and they share one output directory, so an ungated write puts
    one copy of the file per node into the same path -- and since chopper-lib truncates
    rather than appends now, they would overwrite each other's half-written contents
    instead. Counting the header is the cheap way to see either.
    """
    assert _headers(parallel(ranks=ranks, redraw=1), extension) == 1


@requires_mpi
@pytest.mark.parametrize('ranks', MPI_RANKS)
def test_the_grids_carry_the_whole_run(parallel, ranks):
    """The reduced draw count is the run's ncount, not one node's slice of it.

    `counted` takes one unweighted increment per ray that lands in the grid, so with nothing
    absorbed ahead of it the reduced total is exactly NCOUNT however many nodes shared the
    work -- an equality, not a tolerance. Drop the reduction and it becomes NCOUNT / ranks.
    """
    assert _draws(parallel(ranks=ranks, redraw=1)) == pytest.approx(NCOUNT, rel=1e-12)


@requires_mpi
def test_the_intensity_does_not_depend_on_the_rank_count(parallel):
    """Summed, not averaged, and the acceptance left alone.

    The ray weights already carry one over the whole run's ncount, because ESS_butterfly
    reads mcget_ncount() in INITIALIZE and mccode_main slices ncount across the nodes only
    afterwards. So the grids are summed with no division, and the answer is independent of
    how many nodes shared the rays -- within statistics, since the nodes trace different rays.

    A spurious `/= mpi_node_count` shows up here as a factor of `ranks`, and so does an
    acceptance that was reduced along with the grids.
    """
    one = _psd(parallel(ranks=1, redraw=1))
    many = _psd(parallel(ranks=max(MPI_RANKS), redraw=1))

    assert one[1] > 0
    assert many[1] == one[1]                      # the same number of rays reached the monitor
    assert many[0] == pytest.approx(one[0], rel=0.02)
    # and the grids the component writes itself agree the same way
    for extension in ('.total', '.emitted'):
        a = np.loadtxt(parallel(ranks=1, redraw=1).directory / f'guard{extension}', comments='#')
        b = np.loadtxt(parallel(ranks=max(MPI_RANKS), redraw=1).directory / f'guard{extension}',
                       comments='#')
        assert b.sum() == pytest.approx(a.sum(), rel=0.02), extension


@requires_compiler
def test_one_run_writes_one_copy_of_each_file(guard):
    """chopper-lib truncates its output rather than appending to it.

    A caller may write the same grid more than once in a run -- McStas saves on SIGUSR2 and
    carries on, then saves again at the end -- and appending leaves the second copy nose to
    tail with the first, which reads back as one grid of twice the rows. Driving SIGUSR2 from
    a test is awkward; this pins the property the fix rests on, that a write replaces rather
    than extends, which is what makes the repeated save harmless.
    """
    completed = guard(redraw=1)
    for extension in ('.mask', '.total', '.emitted'):
        assert _headers(completed, extension) == 1, extension
