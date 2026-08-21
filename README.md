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

## Describing a chopper with more than one opening

`multi_chopper_parameters` is `{speed, delay, window_count, windows, path}`, where each
of the `windows` is a pair of angles in degrees. An opening edge at angle `a` is on the
beam at `delay + a / (360 * speed)`, and every `1 / |speed|` seconds after that.

Note the sign: `speed` is signed there, and only `|speed|` sets the period. A disk
turning backwards reaches an opening at a positive angle *before* its zero-angle point
rather than after it, so reversing a disk reflects its openings about `delay`. This is
invisible for an opening symmetric about zero -- which is all `single_to_multi_chopper`
builds out of a single-opening chopper -- and matters for every other one.

`multi_chopper_inverse_velocity_windows`, `multi_chopper_inverse_velocity_limits` and
`multi_chopper_wavelength_limits` answer the same questions as their single-opening
counterparts, which they reproduce exactly for a one-window disk.

The mask functions took the sign from `fabs(speed)` before version 3.0.0, so they
mirrored an asymmetric disk that turns backwards. Guard against that where you fill a
`multi_chopper_parameters` whose windows are not symmetric about zero:

```c
#if !defined(CHOPPER_LIB_VERSION) || CHOPPER_LIB_VERSION < 30000
#error "This instrument sets multi-opening choppers; chopper-lib 3.0.0 or newer is required"
#endif
```
