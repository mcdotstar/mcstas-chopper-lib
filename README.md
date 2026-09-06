# A small library for chopper calculations in C, intended for McStas

McStas can be used to simulate time-of-flight neutron scattering instrument,
but was originally develop for constant wavelength instruments.

This small library is intended to provide extra utility to McStas
time-of-flight instruments.

To start, the only utility is intended for use by spectrometers 
which have a train of choppers in their primary spectrometer to define
a wavelength band that can reach the sample.
For direct-geometry spectrometers the band(s) are very sharp,
and for indirect-geometry spectrometers the band(s) are very broad.
In either case, the possible wavelengths that can pass through the
primary spectrometer is reduced by the chopper train.
If the chopper parameters are known, as they must be for the simulation
to progress, they can be used to limit which wavelengths McStas simulates
to only those which *could* make it to the sample position.
The first utility performs the necessary chopper acceptance intersection
calculations and identifies the *envelope* of possible wavelengths.

Embedded in the library code is a utility which returns the *list* of
possible wavelength ranges.
This list could be used as input to a new source which selects
from multiple wavelength bands; or a semi-automatic `GROUP` of 
sources with a pre-source random selection between the bands/sources.
The latter solution forces registration of a new `USERVAR` in 
the particle structure, and a possibly-large group size;
both of which are undesirable.

## Describing a chopper

A chopper is `{speed, delay, angle, path}`: how fast it turns in Hz, when an opening
is on the beam in seconds, how wide that opening is in degrees, and how far it sits
from the source in metres.

`delay` is a time, so it says the same thing whatever the speed is and whichever way
the disk turns — which is how a real chopper is set, and what McStas' `DiskChopper`
acts on. Before version 2.0.0 the second field was a `phase` in degrees, which this
library divided by `360 * fabs(speed)` to recover a delay at every point of use.

Choppers reach the library as flat `double` arrays cast to `chopper_parameters *`, so
that change is invisible to a compiler. Guard against it where you fill the structure:

```c
#if !defined(CHOPPER_LIB_VERSION) || CHOPPER_LIB_VERSION < 20000
#error "This instrument sets chopper delays; chopper-lib 2.0.0 or newer is required"
#endif
```

## Tests

The C library is tested on its own with CTest:

```shell
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

The Python tests under `test/` are separate: they build and run whole McStas
instruments through `niess`, and are driven by pytest.

### Coverage

Configure with `-DCHOPPER_LIB_COVERAGE=ON` and build the `coverage` target. It clears
the counters left by any previous run, runs the suite, and reports which lines of
`chopper-lib.c` the tests reached:

```shell
cmake -S . -B build-coverage -DCHOPPER_LIB_COVERAGE=ON
cmake --build build-coverage --target coverage
```

Instrumentation needs GCC or Clang; configuring with `CHOPPER_LIB_COVERAGE=ON` under
MSVC is an error rather than a silent no-op. The report comes from `gcovr` or `lcov` if
either is installed -- both write a browsable `coverage/index.html`, and `gcovr` also
writes a Cobertura `coverage.xml` for CI -- and otherwise from `gcov`, which ships with
the compiler and gives a per-file summary plus annotated sources under `coverage/gcov`.

## Describing a chopper

`chopper_parameters` is `{speed, delay, beam, edge_count, edges, path}`, and describes a
disk the way the NeXus `NXdisk_chopper` standard and the `CollectorDiskChopper` McStas
component do: `edges` is a flat, increasing list of angles in degrees measured from the
disk's own zero mark, two per opening, and `beam` is the angle of the mark that is on the
beam path at `delay`. An edge at angle `a` is on the beam path at

```
t(a) = delay + (beam - a) / (360 * speed)
```

and every `1 / |speed|` seconds after that. Because only `beam - a` appears, a caller can
hand over the disk's own numbers unchanged: the openings in the frame the disk is drawn
in, and separately where the beam crosses it. There is no conversion to do.

A single-opening chopper of width `w` centred on the mark is `edges = {-w/2, +w/2}` --
negative angles are allowed here, unlike in a NeXus file, because the component allows
them and this library describes the same disk the component does.

Note the sign: `speed` is signed, and only `|speed|` sets the period. A larger angle
reaches the beam *earlier* on a disk turning forwards, so reversing a disk reflects its
openings about `delay`. This is invisible for an opening symmetric about the mark and
matters for every other one.

`chopper_inverse_velocity_windows` lists the openings a train passes;
`chopper_inverse_velocity_limits` and `chopper_wavelength_limits` report the envelope of
that list, with a count so a caller can tell an envelope spanning gaps from a single
window; `chopper_inverse_velocity_time_mask` answers the same question against a
histogram grid.

Version 4.0.0 replaced the `{speed, delay, angle, path}` and
`{speed, delay, window_count, windows, path}` pair of structures with the single one
above, and reversed the sign of the angle term. The field names changed with it, so a
caller written against an older version fails to compile rather than silently placing
every opening on the wrong side of `delay` -- but a caller that fills the structure
positionally does not, and the mask functions took the sign from `fabs(speed)` before
version 3.0.0 besides. Guard where you fill a `chopper_parameters`:

```c
#if !defined(CHOPPER_LIB_VERSION) || CHOPPER_LIB_VERSION < 40000
#error "This instrument describes choppers by edges; chopper-lib 4.0.0 or newer is required"
#endif
```
