from pathlib import Path
from mccode_antlr import Flavor

def compile_and_run(instr,
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

    def actual_compile_run(inside):
        binary, target = mccode_compile(instr, inside, flavor=flavor, **kwargs)
        # The runtime output directory used *can not* exist for McStas/McXtrace to work properly.
        # So find a name inside this directory that doesn't exist (any name should work)
        output = inside / datetime.now().strftime('%Y%m%dT%H%M%S_%f') # isoformat() includes ':' which is not allowed on Windows
        return mccode_run_compiled(binary, target, output, parameters) if run else (None, None)

    if directory is None or not isinstance(directory, Path) or not directory.is_dir():
        with TemporaryDirectory() as tmpdir:
            return actual_compile_run(Path(tmpdir))

    return actual_compile_run(directory)


def make_registry():
    from textwrap import dedent
    from mccode_antlr.reader.registry import InMemoryRegistry
    in_memory_registry = InMemoryRegistry('draft_components')
    comp_name = 'masked_source'
    in_memory_registry.add_comp(comp_name, dedent(r"""
    DEFINE COMPONENT masked_source 
    COPY ESS_butterfly
    DEFINITION PARAMETERS ()
    SETTING PARAMETERS (
      double * choppers, 
      int chopper_count,
      double inverse_velocity_bin,
      double time_bin,
      string filename=0
    )
    OUTPUT PARAMETERS ()
    SHARE %{
      %include "chopper-lib"
    %}
    DECLARE %{
      double * total;
      int * mask;
      unsigned inverse_velocity_count;
      unsigned time_count;
      double minimum_inverse_velocity_edge;
      double inverse_velocity_range;
      double minimum_time_edge;
      double time_range;
    %}
    INITIALIZE %{

      inverse_velocity_range = (Lmax - Lmin) * V2K / 2 / PI;
      if (inverse_velocity_bin <= 0) inverse_velocity_bin = inverse_velocity_range;
      inverse_velocity_count = (unsigned) ceil(inverse_velocity_range / inverse_velocity_bin);
      double * inverse_velocity_edges = (double *) calloc(inverse_velocity_count + 1, sizeof(double));

      time_range = tfocus_width > 0 ? tfocus_width : tmax_multiplier * ESS_SOURCE_DURATION;
      if (time_bin <= 0) time_bin = time_range;
      time_count = (unsigned) ceil(time_range / time_bin);
      double * time_edges = (double *) calloc(time_count + 1, sizeof(double));

      total = (double *) calloc(inverse_velocity_count * time_count, sizeof(double));
      mask = (int *) calloc(inverse_velocity_count * time_count, sizeof(int));

      if (!inverse_velocity_edges || !time_edges || !total || !mask){
        fprintf(stderr, "Out of memory in %s!\n", NAME_CURRENT_COMP);
        exit(1);
      }
      // Use the same mimimum inverse velocity as in ESS_butterfly
      minimum_inverse_velocity_edge = Lmin * V2K / 2 / PI;
      // Use the same minimum time as in ESS_butterfly:
      //    Except we don't know vz yet since we're in initialize, so instead use the
      //    minimum possible vz to give the maximum tfocus_dist/vz, e.g., tfocus_dist*(1/vz)
      //    -> use inverse_velocity_max = minimum_inverse_velocity_edge + inverse_velocity_range:
      // Either all tfocus_* are non-zero or all are zero -- so this is only non-zero if they are all defined:
      minimum_time_edge = tfocus_time - tfocus_dist * (inverse_velocity_range + minimum_inverse_velocity_edge) - tfocus_width/2.0;

      inverse_velocity_edges[0] = minimum_inverse_velocity_edge;
      time_edges[0] = minimum_time_edge;
      for (unsigned iv=0; iv<inverse_velocity_count; ++iv){
        inverse_velocity_edges[iv+1] = inverse_velocity_edges[iv] + inverse_velocity_bin;
        for (unsigned it=0; it<time_count; ++it){
          time_edges[it+1] = time_edges[it] + time_bin;
          mask[iv * time_count + it] = CHOPPER_MASK_INCLUDED;
          total[iv * time_count + it] = 0;
        }
      }

      chopper_parameters * chop_pars = (chopper_parameters *) choppers;

      int allowed_bin_count = chopper_inverse_velocity_time_mask(
        mask, inverse_velocity_count, time_count, 
        inverse_velocity_edges, inverse_velocity_count + 1,
        time_edges, time_count + 1,
        chop_pars, chopper_count,
        1 /* grow the region by one bin in each direction */
      );

      if (allowed_bin_count < 1){
        fprintf(stderr, "Choppers allow no transmission from %s!", NAME_CURRENT_COMP);
        exit(1);
      }
      
      if (!strcmp(filename,"\0")) sprintf(filename,"%s",NAME_CURRENT_COMP);

      if (inverse_velocity_edges) free(inverse_velocity_edges);
      if (time_edges) free(time_edges);
    %}
    TRACE %{
      // Since this is after the trace of ESS_butteryfly, a neutron ray has already been selected.

      // Decide which bin the generated neutron ray should go into:
      double inv_v = 1.0 / sqrt(vx*vx + vy*vy + vz*vz);
      unsigned inverse_velocity_index = (unsigned) floor((inv_v - minimum_inverse_velocity_edge) / inverse_velocity_bin);
      unsigned time_index = (unsigned) floor((t - minimum_time_edge) / time_bin);
      // this indexing must match the internal working of 'chopper_inverse_velocity_time_mask'!
      unsigned linear_index = time_index * inverse_velocity_count + inverse_velocity_index;
      
      printf("Ray at (%f, %f, %f) with v=(%f, %f, %f) -> bin[%d,%d]==%d\n", x, y, z, vx, vy, vz, inverse_velocity_index, time_index, linear_index);

      // Record the total probability distribution sampled
      total[linear_index] += p;

      // If it is masked-out, absorb the ray
      if (mask[linear_index] == CHOPPER_MASK_EXCLUDED) {
        ABSORB;
      } else {
        SCATTER;
      }
    %}
    SAVE %{
      // Now we can find the probability reduction due to our mask
      double relative_probability = chopper_unmasked_probability(total, mask, inverse_velocity_count, time_count);

      int exists = 0;
      FILE * file_mask = mcnew_file(filename, "mask", &exists);
      FILE * file_total = mcnew_file(filename, "total", &exists);

      for (unsigned i=0; i<inverse_velocity_count; ++i){
        for (unsigned j=0; j<time_count; ++j){
          fprintf(file_mask, "%d ", mask[j * inverse_velocity_count + i]);
          fprintf(file_total, "%f ", total[i * time_count + j]);
        }
        fprintf(file_mask, "\n");
        fprintf(file_total, "\n");
      }

//     // Plus save the total probability and mask
//     double * dmask = (double *) calloc(inverse_velocity_count * time_count, sizeof(double));
//     double * zeros = (double *) calloc(inverse_velocity_count * time_count, sizeof(double));
//     for (unsigned idx=0; idx < inverse_velocity_count * time_count; ++idx){
//       dmask[idx] = (double) mask[idx];
//       zeros[idx] = 0.;
//     }
//     DETECTOR_OUT_2D(
//       "Probability Distribution",
//       "Emission Time [s]",
//       "Inverse Velocity [s/m]",
//       minimum_time_edge, minimum_time_edge + time_range,
//       minimum_inverse_velocity_edge, minimum_inverse_velocity_edge + inverse_velocity_range,
//       inverse_velocity_count, time_count,
//       zeros, total, dmask,
//       filename
//     );
//     free(dmask);
//     free(zeros);
    %}
    FINALLY %{
      if (total) free(total);
      if (mask) free(mask);
    %}
    END
    """)
    )

    return comp_name, in_memory_registry


def make_local_registry():
    from mccode_antlr.reader.registry import LocalRegistry
    return LocalRegistry('chopper-lib', '/home/g/Code/mccode/mcstas-chopper-lib')


def insert_source(assembler, source_name):
    from mccode_antlr.common.parameters import InstrumentParameter
    from niess.mccode import ensure_runtime_parameter
    from niess.spatial import mccode_ordered_angles
    from niess.bifrost import Primary
    from niess.components import DiscChopper
    from textwrap import dedent
    from scipp import scalar

    primary = Primary.from_calibration()
    primary.source.n_pulses = 1
    primary.source.accelerator_power = scalar(2.0, unit='MW')

    lines = []
    for x in [y for y, t in primary.items() if t == DiscChopper]:
        comp = getattr(primary, x)
        line = comp.chopper_lib_parameters()
        # hack since we're not building the full instrument, but need the phases and speeds as parameters:
        names = [x.strip() for x in line[1:-1].split(',')]
        for name in names:
            try:
                float(name)
            except ValueError:
                ensure_runtime_parameter(assembler, InstrumentParameter.parse(name))
        lines.append(line + ",")

    assembler.declare(dedent("""
    double * chopper_void;
    unsigned chopper_count;
    """))

    chopper_parameters = "chopper_parameters pars[] ={\n " + '\n '.join(lines) + "\n};"
    void_cast = dedent("""
        chopper_void = (double *) pars;
        chopper_count = sizeof(pars) / sizeof(chopper_parameters);
    """)

    # The BIFROST Primary from_calibration inserts 'source_lambda_min' and 'source_lambda_max'
    # as the wavelength limits of the source component. We can use them and the chopper
    # windows to limit the simulated range in wavelength (-> inverse velocity)
    chopper_calcs = dedent("""
        double lambda_0 = source_lambda_min;
        double lambda_1 = source_lambda_max;
        double pulse_delay = 2.0e-4; // approximate time to peak brightness after protons
        double pulse_width = 2.86e-3; // duration of high flux plateau
        double latest = pulse_delay + pulse_width + 2e-3; // extra for good measure?
        unsigned windows = chopper_wavelength_limits(
         &source_lambda_min, &source_lambda_max, chopper_count, pars, lambda_0, lambda_1, latest
        );
        if (windows == 0){
         printf("Chopper train will not pass wavelengths between %f and %f angstrom\\n", 
                lambda_0, lambda_1);
         printf("Examine the provided chopper speeds and phases.\\n");
        }
        if (windows > 1){
         printf("Chopper train will pass %u wavelength ranges between %f and %f angstrom\\n", 
                windows, lambda_0, lambda_1);
         printf("but only their envelope is considered.\\n");
        }
        printf("Using source lambda limits %f to %f\\n", source_lambda_min, source_lambda_max);
    """)


    assembler.initialize(chopper_parameters + "\n" + void_cast)

    ess_butterfly, parameters = primary.source.__mccode__()
    for f in primary.source.fields():
        if isinstance(p := getattr(primary.source, f), InstrumentParameter):
            ensure_runtime_parameter(assembler, p)

    # We want to use our masked source instead, so update the parameters list:
    parameters['choppers'] = 'chopper_void'
    parameters['chopper_count'] = 'chopper_count'
    parameters['inverse_velocity_bin'] = 0.0001
    parameters['time_bin'] = 0.0001

    at = primary.source.position
    if hasattr(primary.source, 'offset'):
        at += getattr(primary.source, 'offset')
    at = at.to(unit='m').value
    rot = mccode_ordered_angles(primary.source.orientation)

    return assembler.component('source', source_name, at=at, rotate=rot, parameters=parameters)


def make_instr():
    from mccode_antlr.assembler import Assembler
    from mccode_antlr import Flavor
    comp_name, in_memory_registry = make_registry()
    registries = [in_memory_registry, make_local_registry()]
    assembler = Assembler('masked_source_test', flavor=Flavor.MCSTAS, registries=registries)


    origin = assembler.component('origin', 'Progress_bar', at=[0,0,0])
    source = insert_source(assembler, comp_name)

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


def main(e_max, t_psc, source_lambda_min, source_lambda_max, directory):
    instr = make_instr()
    with open(f'{instr.name}.instr', 'w') as file:
        instr.to_file(file)

    args = '-n 1000'
    args += f' source_lambda_min={source_lambda_min} source_lambda_max={source_lambda_max} '
    args += bifrost_chopper_params(e_max, t_psc)
    print(args)
    try:
        compile_and_run(instr, args, run=True, directory=directory)
    except RuntimeError as e:
        print(f'Something went expectedly wrong!\n{e}')


def make_parser():
    from argparse import ArgumentParser

    def is_path(p: str | None):
        if p is None:
            return
        p = Path(p)
        if not p.is_dir():
            p.mkdir(parents=True, exist_ok=False)
        return p

    parser = ArgumentParser()
    parser.add_argument('--e_max', type=int, default=5, help='Maximum incident energy in meV')
    parser.add_argument('--t_psc', type=float, default=0.001, help='PSC opening time in sec')
    parser.add_argument('--lambda_min', type=float, default=0.1, help='Minimum wavelength to simulate')
    parser.add_argument('--lambda_max', type=float, default=10., help='Maximum wavelength to simulate')
    parser.add_argument('--output', type=is_path, default=None, help='Output directory, A temporary directory is used if not provided')
    return parser


if __name__ == '__main__':
    parser = make_parser()
    args = parser.parse_args()
    main(e_max=args.e_max, t_psc=args.t_psc, source_lambda_min=args.lambda_min, source_lambda_max=args.lambda_max, directory=args.output)
    
