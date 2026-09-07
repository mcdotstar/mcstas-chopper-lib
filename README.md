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


## Building the C source

`chopper-lib.c` needs `V2K`, `K2V` and `PI`, and defines none of them. McStas defines
all three in its runtime, and `%include "chopper-lib"` copies the source verbatim into
the generated instrument, so definitions of its own would land there as a second
definition of `PI` beside McStas' -- inert behind include guards, and confusing to
read. Every other build passes them in, and the file refuses to compile without them
rather than quietly fall back to numbers of its own.

The CMake build does that already. A project that pulls this in with FetchContent and
compiles `chopper-lib.c` into a target of its own, rather than linking `chopper_lib`,
passes on the same values:

```cmake
target_compile_definitions(its_target PRIVATE ${CHOPPER_LIB_DEFINITIONS})
```

`CHOPPER_LIB_DEFINITIONS` is cached, so it is readable after `FetchContent_MakeAvailable`
and can be overridden to match a host runtime that defines the constants differently.
Compiling the file by hand takes the same three definitions:

```shell
cc -c chopper-lib.c -DV2K=1.58825361e-3 -DK2V=629.622368 -DPI=3.14159265358979323846
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

`chopper_parameters` is `{speed, delay, beam, edge_count, edges, path, aperture}`, and describes a
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

## A beam of some width

Every field but the last describes a beam of no width: one ray, crossing the disk at the
single angle `beam`. A real beam covers a range of angles, because the opening is angular
and the beam is not, and a neutron crossing `w` metres to one side of the beam centre
reaches an edge `w / d` radians early or late, `d` being the distance from the spindle to
the beam. `aperture` is the width of that range in degrees. Every window widens by half of
it at each end, and nothing else changes: the centre of a window stays where the edges put
it, and a width has no sign, so reversing the disk widens it the same.

It is an angle rather than a width in metres because that is what the disk sees, and
because this library holds no disk geometry to convert one into the other with. The
conversion is not a division either, because a beam window has height as well as width. Its
corners sit further round the disk than its edges do, and the *inner* corners -- nearest the
spindle, where a given width subtends the largest angle -- are furthest of all; while at the
edge of the window the rim has already dropped from `radius` to `sqrt(radius² - (xwidth/2)²)`,
which is where the openings have to hang from if they are to clear the window across its
whole width. For a disk described the way the NeXus standard and McStas' `NXdisk_chopper`
describe one, that leaves the inner corners at `sqrt(radius² - (xwidth/2)²) - yheight` from
the spindle, and

```
aperture = 2 * 180 / pi * atan2(xwidth / 2, sqrt(radius² - (xwidth/2)²) - yheight)
```

Dividing the width by the radius of the beam crossing instead, `xwidth / (radius -
yheight/2)`, misses both corrections and comes out low: 12.7 degrees where the answer is
14.4, for a 100 mm window on a 0.5 m disk with 100 mm openings. `NXdisk_chopper` makes the
same correction to its own geometry when `xwidth` is set, so the two agree about where the
disk is.

Zero -- what a caller that does not set it leaves behind, `aperture` being the last field
-- is the point beam the rest of the fields describe.

This is the field to reach for when a mask cuts a beam a real chopper would pass. The
alternative, growing the finished mask by whole bins, opens the window in inverse velocity
as well as in time, and nothing about a wide beam changes a neutron's wavelength: on one
ESS instrument test with a 100 mm beam on 0.5 m disks, ten bins of growth kept 98.9% of
what the disks passed while opening 27.1% of the frame, where the aperture alone kept 98.0%
of it by opening 23.1% -- and needed no tuning, the number being the geometry.

Note the sign: `speed` is signed, and only `|speed|` sets the period. A larger angle
reaches the beam *earlier* on a disk turning forwards, so reversing a disk reflects its
openings about `delay`. This is invisible for an opening symmetric about the mark and
matters for every other one.

A disk with a `speed` of zero is parked, and is open or shut for good: with no period to
recur on and no delay to apply, all that decides it is whether `beam` falls inside one of
the `edges` pairs. `chopper_parked_is_open` answers that.

A disk parked open constrains nothing -- it has no period, so it passes every inverse
velocity -- and the window and mask functions step over it. A disk parked shut is a beam
stop, and they return nothing at all for one: no windows, no bounds written, no unmasked
bin. That is the answer they already give for a disk with no openings and for a train
whose disks never agree, so a caller has one empty result to handle rather than three
special cases. What an empty result cannot say is which disk emptied it, so a disk parked
shut is named on stdout:

```
chopper-lib: nothing gets through chopper 1, parked with the beam at 90 degrees, where
the disk is solid; the train admits nothing.
```

Scanning a park angle is the case to think about here. At the angles where the disk
blocks the beam there is no band to narrow a source to, and a source that is handed an
empty one has nothing to emit -- which leaves a simulation dividing zero rays by zero
rather than reporting zero intensity from the rays it traced. Decide that at the caller:
either keep the last band that was not empty, or fall back to the full range being
considered, and let the ray tracing put the intensity to zero. Optimising a source
against a disk that blocks the beam is not a thing to do quietly.

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
