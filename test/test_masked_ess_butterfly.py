import time
from contextlib import ContextDecorator
from pathlib import Path
from textwrap import dedent
from pytest import mark
from mccode_antlr import Flavor
from niess.components import ESSource
from scipp import Variable


class Timer(ContextDecorator):
    def __init__(self, name: str = 'block'):
        self.name = name

    def __enter__(self):
        self.start = time.perf_counter()
        return self

    def __exit__(self, *exc):
        self.end = time.perf_counter()
        self.interval = self.end - self.start
        print(f'Timer [{self.name}]: {self.interval:.3f} s')

class MaskedESSource(ESSource):
    identifier_choppers: str
    identifier_chopper_count: str
    inverse_velocity_bin: Variable
    time_bin: Variable

    @classmethod
    def from_source(cls, source: ESSource, cal: dict):
        from scipp import scalar
        source_dict = source.to_dict()
        identifier_choppers = cal.get('identifier_choppers', 'choppers')
        identifier_chopper_count = cal.get('identifier_chopper_count', 'chopper_count')
        inverse_velocity_bin = cal.get('inverse_velocity_bin', scalar(0.001, unit='s/m'))
        time_bin = cal.get('time_bin', scalar(0.001, unit='s'))
        return cls(**source_dict,
                   identifier_choppers=identifier_choppers,
                   identifier_chopper_count=identifier_chopper_count,
                   inverse_velocity_bin=inverse_velocity_bin,
                   time_bin=time_bin)

    def __mccode__(self) -> tuple[str, dict]:
        _, params = super().__mccode__()
        params['choppers'] = self.identifier_choppers
        params['chopper_count'] = self.identifier_chopper_count
        params['inverse_velocity_bin'] = self.inverse_velocity_bin.to(unit='s/m').value
        params['time_bin'] = self.time_bin.to(unit='s').value
        return 'Masked_ESS_butterfly', params


def compiled(method):
    from mccode_antlr.compiler.check import simple_instr_compiles
    if simple_instr_compiles('cc'):
        return method
    @mark.skip(reason=f"Working C compiler required for {method}")
    def skipped_method(*args, **kwargs):
        return method(*args, **kwargs)
    return skipped_method


def time_compile_and_run(instr,
                    parameters,
                    directory: Path | None = None,
                    run=True,
                    dump_source=True,
                    target: dict | None = None,
                    config: dict | None = None,
                    flavor: Flavor = Flavor.MCSTAS):
    from pathlib import Path
    from tempfile import TemporaryDirectory
    from datetime import datetime
    from mccode_antlr.run import mccode_compile, mccode_run_compiled

    kwargs = dict(target=target, config=config, dump_source=dump_source)

    times = {}

    def actual_compile_run(inside):
        with Timer('compile') as t:
            binary, target = mccode_compile(instr, inside, flavor=flavor, **kwargs)
        times['compile_time'] = t.interval
        # The runtime output directory used *can not* exist for McStas/McXtrace to work properly.
        # So find a name inside this directory that doesn't exist (any name should work)
        output = inside / datetime.now().isoformat()
        with Timer('run') as t:
            result = mccode_run_compiled(binary, target, output, parameters) if run else (None, None)
        times['run_time'] = t.interval
        return times, result

    if directory is None or not isinstance(directory, Path) or not directory.is_dir():
        with TemporaryDirectory() as tmpdir:
            return actual_compile_run(Path(tmpdir))

    return actual_compile_run(directory)


def this_registry():
    from git import Repo, InvalidGitRepositoryError
    from mccode_antlr.reader.registry import LocalRegistry
    try:
        repo = Repo('.', search_parent_directories=True)
        root = repo.working_tree_dir
        return LocalRegistry('this_registry', root)
    except InvalidGitRepositoryError as ex:
        raise RuntimeError(f"Unable to identify base repository, {ex}")


def get_registries():
    from mccode_antlr.reader import GitHubRegistry

    registries = [
        'mcstas-frame-tof-monitor', 'mccode-mcpl-filter',
    ]
    registries = [GitHubRegistry(
        name,
        url=f'https://github.com/mcdotstar/{name}',
        filename='pooch-registry.txt',
        version='main'
    ) for name in registries]

    return registries + [this_registry()]


def bifrost_chopper_initialize(assembler, primary, need_pointer: bool = False):
    from niess.components import DiscChopper
    lines = []
    for chopper in [getattr(primary, y) for y, t in primary.items() if t == DiscChopper]:
        line = chopper.chopper_lib_parameters()
        # No need to add phase and speed parameters since we _are_ building the full primary
        lines.append(line)

    assembler.declare(dedent("""
    double * chopper_ptr;
    unsigned chopper_cnt;
    """))
    init = "chopper_parameters pars[] = { " + ", ".join(lines) + " };" + dedent("""
    chopper_ptr = (double *) pars;
    chopper_cnt = sizeof(pars)/sizeof(chopper_parameters);
    
    double lambda_0 = source_lambda_min, lambda_1 = source_lambda_max;
    double pulse_delay = 2.0e-4; // approximate time to peak brightness after protons
    double pulse_width = 2.86e-3; // duration of high flux plateau
    double latest = pulse_delay + pulse_width + 2e-3; // extra for good measure?
    unsigned windows = chopper_wavelength_limits(
     &source_lambda_min, &source_lambda_max, chopper_cnt, pars, lambda_0, lambda_1, latest
    );
    if (windows == 0){
     printf("Chopper train will not pass wavelengths between %f and %f angstrom\\n", lambda_0, lambda_1);
     printf("Examine the provided chopper speeds and phases.\\n");
    }
    if (windows > 1){
     printf("Chopper train will pass %u wavelength ranges between %f and %f angstrom\\n", windows, lambda_0, lambda_1);
     printf("but only their envelope is considered.\\n");
    }
    printf("Using source lambda limits %f to %f\\n", source_lambda_min, source_lambda_max);
    """)
    assembler.initialize(init)


def bifrost_primary(masked: bool = False):
    from scipp import scalar
    from mccode_antlr.assembler import Assembler
    from niess.bifrost import Primary
    name = 'bifrost_masked_ess_source' if masked else 'bifrost_ess_source'
    assembler = Assembler(name, flavor=Flavor.MCSTAS, registries=get_registries())
    primary = Primary.from_calibration()
    primary.source.n_pulses = 1
    primary.source.accelerator_power = scalar(2.0, unit='MW')

    if masked:
        primary.source = MaskedESSource.from_source(primary.source, {
            'identifier_choppers': 'chopper_ptr',
            'identifier_chopper_count': 'chopper_cnt',
            'inverse_velocity_bin': scalar(0.001, unit='s/m'),
            'time_bin': scalar(0.001, unit='s')
        })
    else:
        assembler.initialize('%include "chopper-lib"')

    primary.to_mccode(assembler)

    bifrost_chopper_initialize(assembler, primary, need_pointer=masked)

    return assembler.instrument


def bifrost_chopper_params(e_max, t_psc):
    from chopcal import bifrost
    names = {
        'ps1': 'pulse_shaping_chopper_1',
        'ps2': 'pulse_shaping_chopper_2',
        'fo1': 'frame_overlap_chopper_1',
        'fo2': 'frame_overlap_chopper_2',
        'bw1': 'bandwidth_chopper_1',
        'bw2': 'bandwidth_chopper_2',
    }
    pars = bifrost(e_max, 0, t_psc)
    vals = [f'{v}{y}={getattr(pars[x],y)}' for x, v in names.items() for y in ('phase', 'speed')]
    return ' '.join(vals)


def test_bifrost_ess_source_runs():
    instr = bifrost_primary(masked=False)
    cmds = '-n 1000 source_lambda_min=0.1 source_lambda_max=10 ' + bifrost_chopper_params(5.0, 0.001)
    times, output = time_compile_and_run(instr, cmds, run=True)
    assert times['compile_time'] > 0
    assert times['run_time'] > 0
    assert output is not None


def test_bifrost_masked_ess_source_runs():
    instr = bifrost_primary(masked=True)
    cmds = '-n 1000 source_lambda_min=0.1 source_lambda_max=10 ' + bifrost_chopper_params(5.0, 0.001)
    times, output = time_compile_and_run(instr, cmds, run=True)
    assert times['compile_time'] > 0
    assert times['run_time'] > 0
    assert output is not None


def scipp_monitor_data(output, name='normalization_monitor'):
    import scipp as sc
    dat = output[name]
    x, y, e = dat['t'], dat['I'], dat['I_err']
    rate = sc.array(values=y, variances=e**2, dims=['t'], unit='counts/s')
    return sc.DataArray(rate, coords={'t': sc.array(values= x, dims=['t'], unit='s')})


def test_compare_bifrost_masked_vs_unmasked_ess_source():
    import numpy as np
    instr_masked = bifrost_primary(masked=True)
    instr_unmasked = bifrost_primary(masked=False)
    cmds = '-n 2000000 source_lambda_min=0.1 source_lambda_max=10 ' + bifrost_chopper_params(5.0, 0.001)
    masked_times, (_, output_masked) = time_compile_and_run(instr_masked, cmds, run=True)
    unmasked_times, (_, output_unmasked) = time_compile_and_run(instr_unmasked, cmds, run=True)

    unmasked = scipp_monitor_data(output_unmasked)
    masked = scipp_monitor_data(output_masked)
    diff = (unmasked - masked).sum().data

    from loguru import logger
    mt = masked_times['run_time']
    ut = unmasked_times['run_time']
    logger.info(f"Masked ({mt} s) vs unmasked ({ut} s) difference = {diff:c}")

    assert np.abs(diff.value) < 3*np.sqrt(diff.variance)


def make_parser():
    from argparse import ArgumentParser

    def is_path(pth: str | None):
        if pth is None:
            return None
        pth = Path(pth)
        if not pth.is_dir():
            pth.mkdir(parents=True, exist_ok=False)
        return pth

    p = ArgumentParser()
    p.add_argument('-m','--masked', action='store_true', help='Use masked ESS source')
    p.add_argument('-u','--unmasked', action='store_true', help='Use unmasked ESS source')
    p.add_argument('-n', '--count', type=int, default=1000, help='Number of neutrons to simulate')
    p.add_argument('--e_max', type=int, default=5, help='Maximum incident energy in meV')
    p.add_argument('--t_psc', type=float, default=0.001, help='PSC opening time in sec')
    p.add_argument('--lambda_min', type=float, default=0.1, help='Minimum wavelength to simulate')
    p.add_argument('--lambda_max', type=float, default=10., help='Maximum wavelength to simulate')
    p.add_argument('--output', type=is_path, default=None, help='Output directory, A temporary directory is used if not provided')
    return p


def main(masked, count, e_max, t_psc, source_lambda_min, source_lambda_max, directory):
    cmds = f'-n {count} source_lambda_min={source_lambda_min} source_lambda_max={source_lambda_max} '
    cmds += bifrost_chopper_params(e_max, t_psc)
    print(cmds)
    instr = bifrost_primary(masked=masked)
    time_compile_and_run(instr, cmds, run=True, directory=directory)


if __name__ == '__main__':
    args = make_parser().parse_args()
    main(masked=args.masked and not args.unmasked, count=args.count, e_max=args.e_max, t_psc=args.t_psc,
         source_lambda_min=args.lambda_min, source_lambda_max=args.lambda_max, directory=args.output)

