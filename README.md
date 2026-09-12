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


## The transmitted region

`chopper_inverse_velocity_windows` and `chopper_inverse_velocity_time_mask` each answer a
question about the region a train transmits without building it, and each is wrong in its
own direction. The window function works out every disk's admissible inverse velocities
letting the emission time range over the whole pulse *independently per disk*, then
intersects those ranges — an intersection of projections is a superset of the projection of
the intersection, so it reports bands no single emission time delivers. The mask samples
the region onto a grid, so it misses channels thinner than a bin and counts partly covered
bins whole.

Since 4.2.0 the region itself can be built, exactly, as a set of convex polygons in
(inverse velocity, emission time):

```c
chopper_polygon_set set = chopper_polygon_set_empty();
const chopper_polygon source = chopper_polygon_rectangle(
    inverse_velocity_minimum, inverse_velocity_range, time_minimum, time_range);
chopper_polygon_set_add(&set, &source);
chopper_polygon_set_transmit_train(&set, chopper_count, choppers, NULL);

const double transmitted = chopper_polygon_set_area(&set);   /* s^2/m */
range_set bands = chopper_polygon_set_inverse_velocity_ranges(&set);
```

A neutron emitted at inverse velocity `a` and time `t` reaches path `L` at `t + L*a`, so a
disk open on `[lower, upper]` accepts `lower <= t + L*a <= upper` — a slab between two
parallel lines. A train's acceptance is an intersection of unions of such slabs, and
intersection distributes over union, so the exact acceptance *is* a union of convex pieces,
one per choice of which opening and which turn of each disk a neutron goes through.

Two things follow, and they are what keep this small. Every piece is an intersection of
half-planes, so every piece is convex and the only geometry involved is clipping a convex
polygon by one half-plane — there is no polygon-polygon intersection anywhere. And one clip
adds at most one vertex, so a polygon's vertex count is bounded in advance, which is why
`chopper_polygon` is a fixed-size value that allocates nothing. A rectangle through the six
BIFROST disks comes out as a single five-vertex polygon.

Nothing is removed: the window and mask functions are unchanged and still the right tools
when a grid is what you want, or when the over-estimate is harmless and you would rather
not carry a polygon set around.

`Polygon_ESS_butterfly` is the worked example -- `Masked_ESS_butterfly` over the exact
region rather than a grid of it, and the two are worth running side by side.
`Polygon_ESS_butterfly_image` photographs what each emits. For the two-disc train those
instruments default to, the grid reports an acceptance of 0.290 from 8772 cells where the
region is 0.245 from two triangles: 18% high, and high is the only direction a grid can be
wrong in, since a partly covered cell is kept whole.

What the polygon component gives up is the pictures. `Masked_ESS_butterfly` accumulates
the sampled and emitted probability on its mask grid; there is no grid to accumulate them
on here, and a triangulation is not a sensible thing to histogram into. Two scalar
counters take their place -- rays drawn, and rays landing in the region -- whose ratio is
the acceptance measured from the run rather than from the geometry.

### Sampling it

`chopper_polygon_sampler` does the same job as `chopper_mask_sampler`, in the same shape —
three uniform deviates, one binary search, never rejects — over the exact region rather
than a grid approximation of it:

```c
chopper_polygon_sampler sampler = chopper_polygon_sampler_make(&set, source_area);
chopper_polygon_sampler_draw(&sampler, rand01(), rand01(), rand01(),
                             &inverse_velocity, &time);
p *= sampler.acceptance;
```

`acceptance` is exact here: both areas are known in closed form rather than counted in
cells. The grid sampler can only over-estimate it, because a partly covered cell is
weighted whole — on a BIFROST train over a wide source that is 6.7% high on a
twelve-million-cell grid, and it improves only as fast as the cell count.

### A beam that does not travel in a straight line

A neutron in a guide travels further than the straight line, and how much further depends
on where it bounced. That deviation is in *path*, so what it does to an arrival time is
`deviation * inverse_velocity` — larger for a slow neutron, and nothing at all to the
inverse velocity. So it does not grow the region evenly in every direction; it tilts one of
the two lines bounding each slab, opening it into a wedge.

`chopper_polygon_set_transmit` and `..._transmit_train` take that as a `path_spread` in
metres per disk. Like `aperture` it gives a *support* rather than a distribution: a neutron
is passed if some path in `[path, path + path_spread]` would have got it through, with no
weighting over which. A spread wide enough to reach from one turn of a disk into the next
would stop the transmitted pieces being disjoint and make their areas count twice; that is
refused rather than computed, and a real guide is orders of magnitude inside the limit.

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
instruments, and are driven by pytest. Most need only a C compiler and skip without one;
`test_masked_ess_butterfly.py` also needs `niess` to describe the train, and skips at
collection when it is absent.

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

## Spending a mask rather than throwing rays at it

A mask says which `(inverse velocity, time)` cells a train can pass. The obvious thing to
do with one in a source is to draw a ray, look it up, and absorb it if the cell is
excluded. That is correct -- the discs would have stopped it -- and expensive: a train
passing a few percent of the plane spends the rest of its ray budget on rays that die
where they are born, and the run ends up with the statistics of one a fraction of its
size.

`chopper_mask_sampler` draws from the allowed cells instead. `chopper_mask_sampler_make`
takes a finished mask and the region the caller samples uniformly; `_draw` turns three
uniform deviates into an inverse velocity and a time inside that region, with no rejection
and no loop. The deviates are arguments so the library needs no generator of its own and a
caller inside a McStas TRACE can hand over its own `rand01()`.

```c
chopper_mask_sampler sampler = chopper_mask_sampler_make(
  mask, inverse_velocity_bins, time_bins, inverse_velocity_edges, time_edges,
  inverse_velocity_minimum, inverse_velocity_range, time_minimum, time_range);
...
chopper_mask_sampler_draw(&sampler, rand01(), rand01(), rand01(), &inverse_velocity, &t);
p *= sampler.acceptance;   /* every emitted ray, redrawn or not */
...
chopper_mask_sampler_free(&sampler);
```

That last multiplication is not optional and not a fudge. A ray drawn from a proposal `q`
carrying weight `w` estimates a tally as `E[T] = N E_q[w f]`, and `f` is zero outside the
allowed set because the discs stop those rays. Restricting the draw to the allowed set
multiplies the estimate by `1/Q`, where `Q` is the probability an unrestricted draw lands
there, so multiplying every accepted ray's weight by `Q` puts it back exactly. `acceptance`
is that `Q`. Two numbers it is easy to reach for instead, and neither is right:

- `chopper_unmasked_probability` is the allowed fraction of a *weighted* signal, which is
  the transmission an instrument sees. `Q` counts draws, not intensity.
- a rejection loop's trial count does not yield it either: for a geometric number of trials
  `k`, `E[1/k]` is not `Q`.

`Q` is exact here because it is a ratio of areas, and `_make` clips every cell to the
sampled region before weighting it -- a grid sized with `ceil` runs past that region in its
last row and column, and weighting those cells whole would inflate the answer.

It is exact only while the caller really does draw both coordinates uniformly and
independently of everything else it samples. A source that picks its emission time from a
window centred on the neutron's own velocity does not, and no single factor corrects that.
`Masked_ESS_butterfly` is the worked example: `resample=1` uses the sampler, and its
INITIALIZE refuses time focusing, which is the configuration where the independence fails.

## Writing a grid out under MPI

`chopper_write_mask_to_file` and `chopper_write_total_to_file` truncate what they open, so a
caller that saves the same grid twice in one run -- which McCode does, saving on SIGUSR2 and
carrying on before saving again at the end -- replaces the file rather than leaving the second
copy nose to tail with the first. Before 4.1.0 they appended.

Neither knows anything about MPI, and neither should: a grid is reduced across the nodes by
whoever owns it, and only then written. `Masked_ESS_butterfly` is the worked example again,
and the two halves of it are worth separating because they fail differently.

Every node runs SAVE -- McCode calls `finally()` on all of them and `finally()` calls `save()`
unconditionally -- and `mcuse_dir` hands them all one output directory. So the writing has to
be master's alone, or each file collects one copy per node; and the grids have to be summed
first, or master's copy holds master's share. The sum is a plain `mc_MPI_Sum` with no division
by the node count, because a McCode ray weight already carries one over the *whole* run's
ncount: `mcget_ncount()` returns the full figure in INITIALIZE, and McCode slices it across the
nodes only afterwards.

An `acceptance` is not summed and must not be. It is a function of the mask's geometry, so
every node computes the same number, and reducing it would scale every ray weight in the run by
the node count.

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

Version 4.1.0 added `chopper_mask_sampler` and changed nothing already described, so a
caller that only needs the structures above can keep asking for 4.0.0; one that draws from
a mask should ask for 4.1.0.

```c
#if !defined(CHOPPER_LIB_VERSION) || CHOPPER_LIB_VERSION < 40100
#error "This instrument draws from a chopper mask; chopper-lib 4.1.0 or newer is required"
#endif
```
