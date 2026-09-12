#ifndef CHOPPER_LIB_CHOPPER_LIB_H
#define CHOPPER_LIB_CHOPPER_LIB_H
//
// Created by Gregory Tucker, ESS ERIC on 2023-06-01.
//

/** \file
 *
 * \section versioning Versioning
 *
 * A caller populating a chopper structure is describing a disk in this library's terms,
 * and those terms have changed three times. Before 4.0.0 the structures were flat runs of
 * `double`, so a change to what a field *meant* changed neither size nor layout and an
 * older caller compiled cleanly and computed the wrong answer in silence.
 *
 * `CHOPPER_LIB_VERSION` exists so a caller can refuse to do that. Assert on it wherever
 * a chopper structure is populated:
 *
 *     #if !defined(CHOPPER_LIB_VERSION) || CHOPPER_LIB_VERSION < 40000
 *     #error "This instrument describes disks by their slit edges; chopper-lib 4.0.0 or newer is required"
 *     #endif
 *
 * From 4.0.0 the structure holds a pointer and a count, so most mistakes are now compile
 * errors rather than silent ones -- but the guard still earns its place, because the
 * *meaning* of an edge angle changed at the same time and no compiler can see that.
 *
 * The major version changes when the meaning or layout of a structure changes.
 *
 * 4.2.1
 *     `range_set_sort` merges correctly. It lost the extent of a range containing the one
 *     after it, and its answer depended on the order `qsort` left tied lower edges in --
 *     which is not fixed across platforms, so a train's admitted band could differ
 *     between Windows and Linux. Present since 2.0.0 and reached whenever the ranges
 *     overlap enough to tie, which a disk `aperture` makes likely.
 *
 *     `range_sort` is gone. It took its `range` by value, so it swapped the edges of a
 *     copy and returned -- it never did anything, whatever its name and its documentation
 *     said. Nothing called it: `range_set_sort` puts its own ranges the right way round.
 *     A removal rather than a fix because the signature is the fault, and a removal in a
 *     patch release because the rule above is about what a structure means, and a
 *     function that did nothing cannot have meant anything.
 *
 * 4.2.0
 *     `chopper_polygon` and the functions around it carry the transmitted region of
 *     (inverse velocity, time) as a set of convex polygons, which is what a chopper train
 *     actually passes rather than an approximation of it. `chopper_inverse_velocity_windows`
 *     and `chopper_inverse_velocity_time_mask` are unchanged and still here; see the note
 *     on `chopper_polygon_set_transmit_train` for what each of them gets wrong and why.
 *     Additive: nothing already described means anything different.
 *
 * 4.1.0
 *     `chopper_mask_sampler` draws an (inverse velocity, time) pair from the allowed
 *     cells of a finished mask, and carries the `acceptance` a caller has to multiply
 *     the ray weight by to keep the sampling unbiased. Additive: nothing already
 *     described means anything different.
 *
 * 4.0.0
 *     One `chopper_parameters` describes a disk of any number of openings, by its slit
 *     edges. `chopper_window`, `multi_chopper_parameters`, `single_to_multi_chopper` and
 *     the parallel `multi_chopper_*` functions are gone; the `angle` field is gone with
 *     them.
 *
 *     `edges` is the `slit_edges` of the NeXus NXdisk_chopper specification: angles from
 *     the disk's top-dead-centre mark, strictly increasing, opening edge first. Window
 *     angles were measured against the beam and in the opposite sense, so an edge that
 *     was `a` is now `beam - a` and each pair has swapped ends.
 *
 *     The new `beam` field is the angle from the mark to where the beam crosses the disk
 *     -- the `beam` of the same NeXus specification. It used to have to be folded into
 *     `delay` or into every window angle by the caller; it is a field now, so a disk is
 *     described the same way here as it is in the file and in the McStas component.
 *
 *     `chopper_parameters` ends with an `aperture`: the angular width of the beam on the
 *     disk, in degrees, which widens every window it computes by half of it at each end.
 *     It is the last field, so a caller filling the structure positionally without one
 *     still compiles and still describes the point beam it described before.
 *
 *     A disk parked shut is a beam stop. Every disk with a `speed` of zero used to be
 *     left out of the calculation, which is right for one parked open -- it has no
 *     period, so it constrains nothing -- and wrong for one parked shut, which passes
 *     nothing at any time and was reported as a band the instrument would not deliver.
 *     The window and mask functions now empty their answer for one, as they already did
 *     for a disk with no openings, and name it on stdout on the way past.
 *     `chopper_parked_is_open` is the predicate they use, and says whether a disk that is
 *     not turning stands open on the beam.
 * 3.0.0
 *     A `multi_chopper_parameters` window angle is placed with the *signed* `speed`:
 *     an opening at angle `a` is on the beam at `delay + a / (360 * speed)`. The mask
 *     functions previously used `fabs(speed)`, which reflects an asymmetric disk about
 *     its `delay` when the disk turns backwards. Only a `windows` array symmetric about
 *     zero -- all `single_to_multi_chopper` produces -- is unaffected, so a caller that
 *     compensated for the old behaviour now has the error twice over.
 * 2.0.0
 *     `chopper_parameters` and `multi_chopper_parameters` take a `delay` in seconds
 *     where they previously took a `phase` in degrees. See the note on `delay` below.
 * 1.0.0
 *     Unversioned releases, taking `phase`.
 */
#define CHOPPER_LIB_VERSION_MAJOR 4
#define CHOPPER_LIB_VERSION_MINOR 2
#define CHOPPER_LIB_VERSION_PATCH 1
/** Single integer form, MAJOR*10000 + MINOR*100 + PATCH, for comparison in `#if` */
#define CHOPPER_LIB_VERSION (CHOPPER_LIB_VERSION_MAJOR * 10000 \
                           + CHOPPER_LIB_VERSION_MINOR * 100 \
                           + CHOPPER_LIB_VERSION_PATCH)

/** A contiguous range characterized by two edge values
 *
 * @param minimum The lower edge of the contiguous range
 * @param maximum The upper edge of the contiguous range
 */
struct range_struct {
  double minimum;
  double maximum;
};
typedef struct range_struct range;
struct int_range_struct {
  int minimum;
  int maximum;
};
typedef struct int_range_struct int_range;

/** Determine whether two ranges overlap, and if so, characterize their overlapping type
 *
 * @param a A first range
 * @param b A second range
 * @return 0 if the two ranges do not overlap, 1 if they are identical, +-2 if the first range extends to higher values,
 *         +-3 if the second range extends to higher values, (+) if one range is inside the other, and (-) if the
 *         overlapping region is a subset of both the first and second ranges.
 */
int classify_range_overlap(const range * a, const range * b);
/** Compare two ranges based only on their minimum edges
 *
 * @param a A pointer to the first range
 * @param b A pointer to the second range
 * @return 1 if the lower edge of the first range is higher than that of the second range, -1 if the second range
 *         lower edge is higher, or 0 if the two lower edges are equivalent.
 */
int compare_ranges(const range * a, const range * b);

/** A collection of contiguous regions which are not fully contiguous */
struct range_set_struct {
  unsigned count;
  range * ranges;
};
typedef struct range_set_struct range_set;

/** Sort the contiguous sub-ranges of a set of ranges
 *
 * @param s The set of ranges to sort
 * @return The sub-ranges sorted by lower edge with overlapping or contiguous regions merged
 * @warning The function *may* allocate a new `range_set` for output or return the input. In either case the input
 *          structure `.ranges` field is sorted. Care should be taken to release memory from the returned value.
 */
range_set range_set_sort(range_set s);
/** Find the intersection of two sets of contiguous ranges
 *
 * @param ain The first set of ranges
 * @param bin The second set of ranges
 * @return A set of ranges covered by both input range sets. A range is output if and only if it is part of *both*
 *         of the input range sets.
 */
range_set range_intersection(range_set ain, range_set bin);


/** The parameters of a disk chopper, of one opening or several.
 *
 * @param speed The rotation speed of the disk in Hz; negative turns it the other way
 * @param delay When the disk point at angle `beam` is on the path, in seconds
 * @param beam The angle from the disk's zero mark to where the beam crosses it, in degrees
 * @param edge_count The number of entries in `edges`: two per opening, so always even
 * @param edges The opening and closing edge of each opening, in degrees from the zero mark
 * @param path The path length from the 'zero'-time source to the disk position, in meters
 * @param aperture The angular width of the beam where it crosses the disk, in degrees
 *
 * `edges` is the `slit_edges` of the NeXus NXdisk_chopper specification, and what McStas'
 * `CollectorDiskChopper` takes: an even number of angles measured from the top-dead-centre
 * mark, strictly increasing, the opening edge of each slit first, spanning less than one
 * turn. A slit straddling the mark is written with a final edge past 360 -- `{350, 370}`
 * rather than `{350, 10}` -- so the pairs stay ordered and each width is a difference.
 *
 * Nothing here requires the angles to be positive. NXdisk_chopper does, and a caller
 * writing a NeXus file should rotate which slit comes first to satisfy it, but a disk with
 * one opening astride its mark reads better here as `{-85, 85}` than as `{275, 445}`.
 *
 * An edge at angle `a` is on the beam at
 *
 *     t(a) = delay + (beam - a) / (360 * speed)
 *
 * and every `1 / |speed|` seconds thereafter. Two things follow from the signs. The
 * angular term is negated, because `edges` increase in the direction the NeXus
 * specification measures and a disk carries a larger angle *towards* the beam; and it
 * keeps the sign of `speed`, so reversing the disk brings an opening onto the beam on the
 * other side of `delay`.
 *
 * `delay` is a time, so it is unaffected by the sign of `speed`, and openings recur at
 * `delay + n / speed` for integer `n`. It is what McStas' `DiskChopper` acts on and what a
 * real chopper is set with. It was a `phase` in degrees before version 2.0.0, from which
 * this library recovered a delay by dividing by `360 * fabs(speed)` at every point of use;
 * a delay says the same thing without needing to know the speed, and -- unlike a phase,
 * which wraps at one revolution -- may exceed a single period.
 *
 * A disk of `n` identical, evenly spaced openings may be described as one opening turning
 * `n` times as fast, if that is more convenient; `delay` is unaffected by the substitution.
 *
 * Every other field describes a beam of no width: one ray, crossing the disk at the single
 * angle `beam`. A real beam covers a range of angles, because an opening is angular and the
 * beam is not -- a neutron that crosses `w` metres to one side of the beam centre reaches
 * the same opening edge `w / d` radians early or late, where `d` is the distance from the
 * spindle to the beam. `aperture` is the width of that range in degrees, so an opening is
 * open half of it either side of the times the edges alone would give. Zero, the value a
 * caller that does not set it leaves behind, is the point beam the rest of the fields
 * describe.
 *
 * It is an angle rather than a width in metres because that is what the disk sees, and
 * because this library holds no disk geometry to convert one into the other with. The
 * conversion is not a division either, because a beam window has height as well as width:
 * its corners are further round the disk than its edges are, and its *inner* corners --
 * nearest the spindle, where a given width is the largest angle -- are the furthest of
 * all. For a disk of `radius` whose openings clear a window `xwidth` wide by `yheight`
 * radially, described the way the NeXus standard and McStas' NXdisk_chopper describe one,
 * the rim has dropped to `sqrt(radius^2 - (xwidth/2)^2)` at the edge of the window, so the
 * openings reach in to that less `yheight`, and
 *
 *     aperture = 2 * 180 / pi * atan2(xwidth / 2, sqrt(radius^2 - (xwidth/2)^2) - yheight)
 *
 * Taking the width over the radius of the beam crossing instead, `xwidth / (radius -
 * yheight/2)`, misses both corrections and comes out low -- 12.7 degrees where the real
 * figure is 14.4, for a 100 mm window on a 0.5 m disk with 100 mm openings.
 *
 * A wider aperture opens each window in time, and only in time -- nothing about a wide beam
 * changes a neutron's wavelength -- which is what makes it the right place to account for
 * the width. Growing a finished mask by whole bins does the same job in both directions at
 * once, and admits inverse velocities no disk ever passes.
 */
struct chopper_parameters_struct {
  double speed; // rotation frequency in Hz
  double delay; // when the disk point at angle `beam` is on the path, in seconds
  double beam; // from the zero mark to the beam crossing, in degrees
  unsigned edge_count; // number of entries in edges, two per opening
  double * edges; // opening and closing edge of each opening, in degrees from the mark
  double path; // average(?) path length from source to this chopper in meters
  double aperture; // angular width of the beam on the disk, in degrees; 0 is a point
};
typedef struct chopper_parameters_struct chopper_parameters;


/** Whether a disk that is not turning stands open on the beam
 *
 * A parked disk is open or shut for good: with no speed there is no period to recur on
 * and no delay to apply, so either `beam` is inside one of the `edges` pairs or it is on
 * the solid part of the disk. Angles fold, so an opening written across the mark --
 * `{350, 370}` -- and a negative `beam` both work.
 *
 * The window and the mask functions leave a disk parked open out of their calculation,
 * because a disk with no period constrains no inverse velocity. One parked shut is a
 * beam stop: they name it on stdout and return nothing at all, the same answer they give
 * for a disk with no openings and for a train whose disks never agree. Call this first
 * to tell that apart from a train that is merely over-constrained, or to decide what a
 * shut disk should mean in your own terms.
 *
 * @param chopper The disk to test; its `speed` is not read, since a turning disk stands
 *                open on the beam once a period whatever its angles are
 * @return 1 if the beam crosses an opening, 0 if it crosses the disk body or the disk has
 *         no openings at all
 */
int chopper_parked_is_open(chopper_parameters chopper);

/** Find the possible inverse velocity window(s) that are admitted by a series of disk choppers
 *
 * @param count The number of disk choppers provided
 * @param choppers The parameters of the disk choppers
 * @param inv_v_min The minimum inverse velocity to be considered -- likely matching a guide cutoff
 * @param inv_v_max The maximum inverse velocity to be considered -- how long before a neutron is no-longer interesting
 * @param latest_emission How long after time-zero can a neutron start its journey, effects minimum inverse velocities
 * @return One or more inverse velocity ranges that can pass through the chopper train as a `range_set`
 * @note A disk parked open is left out: it has no period, so it constrains no inverse
 *       velocity. One parked shut passes nothing at any time, so the returned set is
 *       empty and the disk is named on stdout -- see `chopper_parked_is_open`.
 * @warning The returned value's `ranges` property is allocated in the function and must be freed at calling scope.
 */
range_set chopper_inverse_velocity_windows(unsigned count, const chopper_parameters * choppers,
                                           double inv_v_min, double inv_v_max, double latest_emission);

/** Find the enveloping limits of the possible inverse velocity window(s) that are admitted by a chopper train
 *
 * @param lower Output lower inverse velocity limit, only set if the return value is finite
 * @param upper Output upper inverse velocity limit, only set if the return value is finite
 * @param count The number of choppers in the train
 * @param choppers Parameters for each chopper
 * @param inv_v_min The minimum inverse velocity to be considered
 * @param inv_v_max The maximum inverse velocity to be considered
 * @param latest_emission How long after time-zero a neutron can start along the flight path
 * @return The number of inverse velocity windows admitted by the choppers, if greater than one the lower and upper
 *         values include in their range inverse velocities which are not passed by the chopper train.
 */
unsigned chopper_inverse_velocity_limits(double * lower, double * upper,
                                         unsigned count, const chopper_parameters * choppers,
                                         double inv_v_min, double inv_v_max, double latest_emission);

/** Find the enveloping limits of the possible wavelength window(s) that are admitted by a chopper train
 *
 * @param lower Output lower wavelength limit, only set if the return value is finite
 * @param upper Output upper wavelength limit, only set if the return value is finite
 * @param count The number of choppers in the train
 * @param choppers Parameters for each chopper
 * @param lambda_min The minimum wavelength to be considered
 * @param lambda_max The maximum inverse velocity to be considered
 * @param latest_emission How long after time-zero a neutron can start along the flight path
 * @return The number of windows admitted by the choppers, if greater than one the lower and upper
 *         values include in their range wavelengths which are not passed by the chopper train.
 */
unsigned chopper_wavelength_limits(double * lower, double * upper, unsigned count, const chopper_parameters * choppers,
                                   double lambda_min, double lambda_max, double latest_emission);

/** Create a mask of allowed (inverse_velocity, time) bins based on chopper parameters
 *
 * @param mask [out] An array of integers to be filled with 1 (allowed) or 0 (blocked), of size (inverse_velocity_count-1) * (time_count-1)
 * @param mask_inverse_velocity_count [in] The number of inverse velocity bins in the mask (should be inverse_velocity_count - 1)
 * @param mask_time_count [in] The number of time bins in the mask (should be time_count - 1)
 * @param inverse_velocities [in] An array of inverse velocities (in s/m), of size inverse_velocity_count
 * @param inverse_velocity_count [in] The number of inverse velocities provided
 * @param times [in] An array of times at the source position (in s), of size time_count
 * @param time_count [in] The number of times provided
 * @param choppers [in] An array of chopper parameters, of size chopper_count
 * @param chopper_count [in] The number of choppers provided
 * @param grow_mask [in] Expand the allowed regions by this number of bins in each direction
 * @return The number of unmasked (allowed) (inverse_velocity, time) bins
 * @note A disk parked open is left out, as it is by the window functions; one parked
 *       shut masks off every bin and is named on stdout -- see `chopper_parked_is_open`.
 */
unsigned chopper_inverse_velocity_time_mask(
  int * mask, unsigned mask_inverse_velocity_count, unsigned mask_time_count,
  const double * inverse_velocities, unsigned inverse_velocity_count,
  const double * times, unsigned time_count,
  const chopper_parameters * choppers, unsigned chopper_count,
  int grow_mask
  );

/** \brief Calculate the relative probability in the unmasked regions of a signal
 *
 * @param signal The probability to consider, as a flattened 2D array of size mask_inverse_velocity_count * mask_time_count
 * @param mask A mask of allowed (1) and disallowed (0) bins, of size mask_inverse_velocity_count * mask_time_count
 * @param mask_inverse_velocity_count The number of inverse velocity bins in the mask
 * @param mask_time_count The number of time bins in the mask
 * @return The mask expectation value of the signal, i.e., the sum of the signal in allowed bins divided by the total signal
 */
double chopper_unmasked_probability(
  const double * signal, const int * mask, unsigned mask_inverse_velocity_count, unsigned mask_time_count
  );

enum mask_values {
  CHOPPER_MASK_EXCLUDED = 0,
  CHOPPER_MASK_INCLUDED = 1,
  CHOPPER_MASK_GROWN = 100
};

/** A direct sampler over the allowed cells of a finished mask
 *
 * A source that draws an (inverse velocity, time) pair and throws it away when the mask
 * excludes it spends its whole ray budget to keep the fraction the choppers pass. Drawing
 * from the allowed cells in the first place keeps all of it, and is the same distribution
 * -- provided the weight is corrected, which is what `acceptance` is for.
 *
 * The correction is not a matter of taste. A ray drawn from proposal `q` and carrying
 * weight `w` estimates a downstream tally as `E[T] = N E_q[w f]`, and `f` is zero outside
 * the allowed set `A` because the discs stop those rays. Drawing from `q` restricted to
 * `A` instead multiplies that by `1 / Q`, where
 *
 *     Q = P(a draw from q lands in A)
 *
 * so multiplying every accepted ray's weight by `Q` puts it back. Two things this `Q` is
 * not, both of which are easy to reach for:
 *
 *   - It is not `chopper_unmasked_probability`. That is the allowed fraction of a
 *     *weighted* signal, which is the transmission an instrument sees and the wrong
 *     normalisation for this. `Q` counts draws, not intensity.
 *   - It is not recoverable from a rejection loop's trial count. For a geometric number
 *     of trials `k`, `E[1/k]` is not `Q`, so per-ray attempt counting biases the answer.
 *
 * `Q` is exact here because it is a ratio of areas: `chopper_mask_sampler_make` is told
 * the region the caller samples uniformly, clips every cell to it, and returns the
 * allowed area over the whole area. It is exact only while the caller really does sample
 * both coordinates uniformly and independently of everything else it draws -- see the
 * note on `chopper_mask_sampler_make`.
 *
 * Drawing costs one binary search and three uniform deviates, and never rejects.
 *
 * @param count The number of allowed cells the sampler draws from
 * @param acceptance `Q`: the allowed area over the sampled area, and the factor a ray
 *                   weight must be multiplied by
 * @param cumulative `count` entries increasing to 1, the area-weighted cell distribution
 * @param inverse_velocity_low The low edge of each cell, clipped to the sampled region
 * @param inverse_velocity_width Its width after clipping; never negative, never zero
 * @param time_low The low time edge of each cell, clipped the same way
 * @param time_width Its width after clipping
 */
struct chopper_mask_sampler_struct {
  unsigned count;
  double acceptance;
  double * cumulative;
  double * inverse_velocity_low;
  double * inverse_velocity_width;
  double * time_low;
  double * time_width;
};
typedef struct chopper_mask_sampler_struct chopper_mask_sampler;

/** Build a sampler over the allowed cells of a mask
 *
 * The mask grid and the region a caller samples are not the same rectangle. A grid sized
 * with `ceil` runs past the region in its last row and column, and a cell there is only
 * partly reachable; weighting it whole would over-represent it and inflate `acceptance`.
 * Every cell is therefore clipped to
 * `[inverse_velocity_minimum, inverse_velocity_minimum + inverse_velocity_range]` and
 * `[time_minimum, time_minimum + time_range]` before it is weighted, and a cell left with
 * no area is dropped.
 *
 * @param mask The finished mask, as `chopper_inverse_velocity_time_mask` leaves it:
 *             `CHOPPER_MASK_EXCLUDED` where nothing passes and `CHOPPER_MASK_INCLUDED`
 *             everywhere else, grown cells included
 * @param mask_inverse_velocity_count Inverse velocity bins in the mask
 * @param mask_time_count Time bins in the mask
 * @param inverse_velocities The `mask_inverse_velocity_count + 1` bin edges, increasing
 * @param times The `mask_time_count + 1` time bin edges, increasing
 * @param inverse_velocity_minimum The low edge of the region the caller samples
 * @param inverse_velocity_range Its width; the caller draws uniformly across it
 * @param time_minimum The low edge of the time region the caller samples
 * @param time_range Its width; the caller draws uniformly across it too
 * @return A sampler whose `cumulative` and edge arrays are allocated here and must be
 *         released with `chopper_mask_sampler_free`. A mask that allows nothing, or
 *         allows nothing inside the sampled region, comes back with `count` 0,
 *         `acceptance` 0 and no allocations; drawing from it is a caller error.
 *
 * @note The `acceptance` this returns is only the right weight correction while the two
 *       coordinates are drawn uniformly, independently of each other, and independently
 *       of everything else the caller samples. A source that picks its emission time from
 *       a window centred on the neutron's own velocity -- McStas' `ESS_butterfly` under
 *       time focusing does exactly that -- breaks the independence, and no single number
 *       corrects it.
 */
chopper_mask_sampler chopper_mask_sampler_make(
  const int * mask, unsigned mask_inverse_velocity_count, unsigned mask_time_count,
  const double * inverse_velocities, const double * times,
  double inverse_velocity_minimum, double inverse_velocity_range,
  double time_minimum, double time_range
  );

/** Zero a sampler, so it can be freed or tested before it has been built
 *
 * A caller that builds a sampler only on some paths still has to be able to free it on all
 * of them. This puts one in the state `chopper_mask_sampler_free` leaves behind: `count`
 * and `acceptance` zero, every pointer NULL.
 *
 * @param sampler The sampler to empty; nothing it currently points at is released, so do
 *                not call this on a built sampler in place of `chopper_mask_sampler_free`
 */
void chopper_mask_sampler_empty(chopper_mask_sampler * sampler);

/** Release what `chopper_mask_sampler_make` allocated, and leave an empty sampler behind
 *
 * @param sampler The sampler to empty; a NULL pointer, or one already emptied, is fine
 */
void chopper_mask_sampler_free(chopper_mask_sampler * sampler);

/** Draw one (inverse velocity, time) pair uniformly from the allowed region
 *
 * The three deviates are arguments rather than drawn here so that the library needs no
 * random number generator of its own, and so that a caller inside a McStas TRACE can
 * hand over `rand01()` and stay on whatever generator the instrument was built with.
 *
 * @param sampler A sampler with a non-zero `count`
 * @param cell_deviate A uniform deviate on [0, 1), choosing which allowed cell
 * @param inverse_velocity_deviate A uniform deviate on [0, 1), placing the point across it
 * @param time_deviate A uniform deviate on [0, 1), placing the point up it
 * @param inverse_velocity [out] The drawn inverse velocity, in s/m
 * @param time [out] The drawn time, in s
 *
 * @note The caller still owes the ray weight a factor of `sampler->acceptance`.
 */
void chopper_mask_sampler_draw(
  const chopper_mask_sampler * sampler,
  double cell_deviate, double inverse_velocity_deviate, double time_deviate,
  double * inverse_velocity, double * time
  );

/*************************** transmitted phase space *********************************/

/** \section polygons The transmitted region, exactly
 *
 * `chopper_inverse_velocity_windows` and `chopper_inverse_velocity_time_mask` each answer
 * a question about the transmitted region without ever building it, and each is wrong in
 * its own direction. The window function works out every disk's admissible inverse
 * velocities letting the emission time range over the whole pulse *independently per
 * disk*, then intersects those ranges; an intersection of projections is a superset of the
 * projection of the intersection, so it reports bands that no single emission time
 * delivers. The mask samples the region onto a grid, so it both misses channels thinner
 * than a bin and counts partly covered bins whole.
 *
 * The region itself is not hard to build. A neutron emitted at inverse velocity `a` and
 * time `t` reaches path `L` at `t + L*a`, so a disk open on `[lower, upper]` accepts
 *
 *     lower <= t + L*a <= upper
 *
 * -- a slab between two parallel lines. A train's acceptance is an intersection of unions
 * of such slabs, one union per disk, and intersection distributes over union, so the exact
 * acceptance *is* a union of convex pieces: one per choice of which opening and which turn
 * of each disk a neutron goes through. Nothing here approximates anything.
 *
 * Two consequences make it much less code than it sounds. Every piece is an intersection
 * of half-planes, so every piece is convex and the only geometry needed is clipping a
 * convex polygon by one half-plane -- there is no polygon-polygon intersection anywhere in
 * this file. And a clip adds at most one vertex, so a polygon's vertex count is bounded in
 * advance, which is why `chopper_polygon` can be a fixed-size value with no allocation of
 * its own.
 *
 * The pieces are disjoint whenever a disk is shut for part of every turn, so
 * `chopper_polygon_set_area` is a plain sum and a sampler may pick a piece by area without
 * double counting. `chopper_polygon_set_transmit` checks the one case that would break
 * that; see its note on `path_spread`.
 */

/** The most vertices one polygon can reach.
 *
 * A clip against one half-plane adds at most one vertex and a disk costs two clips, so a
 * `V`-vertex source polygon through `n` disks needs `V + 2n`. 32 covers a rectangle
 * through fourteen disks, and a six-disk train from a rectangle has been measured to stay
 * at six. Exceeding it is an error rather than a truncation: dropping a vertex would
 * quietly *enlarge* the region, which is the class of mistake this whole section exists to
 * end.
 */
#define CHOPPER_POLYGON_MAX_VERTICES 32

/** How thin a polygon may be before it is treated as nothing.
 *
 * Relative to the polygon's own bounding box, so it is free of any choice of units --
 * which matters more than it looks. The same region is 1.3e-07 in s^2/m and 1.3e+02 in
 * ms^2/km, so an absolute tolerance means different things in different unit bases, and a
 * library that took one would give different answers to callers working in different
 * units. This one does not.
 *
 * Clipping along an edge a polygon already has leaves a piece of no area, and those are
 * what this is for. A genuinely narrow transmission channel is nowhere near it: a sliver
 * one part in a million of its own length still has an area a millionth of its bounding
 * box, six orders above this.
 */
#define CHOPPER_POLYGON_AREA_TOLERANCE 1e-12

/** A point of the (inverse velocity, time) plane: s/m and s, at the source. */
struct chopper_point_struct {
  double inverse_velocity;
  double time;
};
typedef struct chopper_point_struct chopper_point;

/** A convex region of (inverse velocity, emission time), its vertices in order.
 *
 * Fixed capacity, so this is a value: copy it, put it on the stack, let it go out of
 * scope. It owns nothing and needs no freeing. Only `chopper_polygon_set` allocates.
 */
struct chopper_polygon_struct {
  unsigned count;
  chopper_point vertex[CHOPPER_POLYGON_MAX_VERTICES];
};
typedef struct chopper_polygon_struct chopper_polygon;

/** A union of convex regions -- what a source emits, or what a train transmits.
 *
 * Disjoint as `chopper_polygon_set_transmit` leaves it, because two turns of one disk
 * cannot both pass the same neutron. That is what lets `chopper_polygon_set_area` add
 * rather than needing inclusion-exclusion.
 */
struct chopper_polygon_set_struct {
  unsigned count;
  unsigned capacity;
  chopper_polygon * polygon;
};
typedef struct chopper_polygon_set_struct chopper_polygon_set;

/** The rectangle a source emits into.
 *
 * @param inverse_velocity_minimum Low edge, s/m
 * @param inverse_velocity_range Its width; must be positive
 * @param time_minimum Low edge of the emission window, s
 * @param time_range Its width; must be positive
 * @return The rectangle, counter-clockwise; a polygon of no vertices if either range is
 *         not positive
 */
chopper_polygon chopper_polygon_rectangle(double inverse_velocity_minimum,
                                          double inverse_velocity_range,
                                          double time_minimum, double time_range);

/** Area by the shoelace formula, in s^2/m. Zero for fewer than three vertices. */
double chopper_polygon_area(const chopper_polygon * polygon);

/** The range of `alpha*inverse_velocity + beta*time` over a polygon.
 *
 * With `alpha` a path and `beta` 1 this is the earliest and latest a neutron in the
 * polygon reaches that path, which is how the turns of a disk that could matter are found.
 *
 * @param polygon The polygon to measure; may be empty, in which case nothing is written
 * @param alpha Weight on inverse velocity
 * @param beta Weight on time
 * @param lower [out] The smallest value, if not NULL
 * @param upper [out] The largest, if not NULL
 */
void chopper_polygon_extent(const chopper_polygon * polygon, double alpha, double beta,
                            double * lower, double * upper);

/** Keep the part of a polygon where `alpha*inverse_velocity + beta*time <= c`.
 *
 * The only geometry in this file; everything above and below is built from it. A disk at
 * path `L` is `alpha = L, beta = 1`; an emission window is `alpha = 0`; a bound on inverse
 * velocity alone is `beta = 0`.
 *
 * @param polygon The polygon, clipped in place. Left with no vertices if nothing survives.
 * @param alpha Weight on inverse velocity
 * @param beta Weight on time
 * @param c The bound
 * @return 0 if the result would need more than `CHOPPER_POLYGON_MAX_VERTICES` vertices,
 *         leaving `polygon` untouched; 1 otherwise
 */
int chopper_polygon_clip_halfplane(chopper_polygon * polygon,
                                   double alpha, double beta, double c);

/** Keep what a disk passes when the path to it is known only to lie in a range.
 *
 * A neutron in a guide travels further than the straight line between two points, and how
 * much further depends on where it bounced. That deviation is in *path*, so what it does
 * to an arrival time is `deviation * inverse_velocity` -- larger for a slow neutron, and
 * nothing at all to the inverse velocity itself. So it does not grow the polygon evenly in
 * every direction; it tilts one of the two bounding lines.
 *
 * A neutron emitted at `(a, t)` arrives somewhere in `[t + shortest*a, t + longest*a]`, and
 * passes if any of that lands in the window:
 *
 *     t + shortest_path*a <= upper     and     t + longest_path*a >= lower
 *
 * Two half-planes again, of different slopes, so the slab opens into a wedge as the
 * inverse velocity grows. Equal paths give the ordinary slab.
 *
 * @param polygon The polygon, clipped in place
 * @param shortest_path The shortest path a neutron could have taken to this disk, m
 * @param longest_path The longest; must not be less than `shortest_path`
 * @param beta Weight on time, normally 1
 * @param lower The window opens
 * @param upper The window closes
 * @return 0 on vertex overflow, 1 otherwise
 */
int chopper_polygon_clip_wedge(chopper_polygon * polygon, double shortest_path,
                               double longest_path, double beta,
                               double lower, double upper);

/** An empty set, owning nothing. Safe to free, and what a failed call leaves behind. */
chopper_polygon_set chopper_polygon_set_empty(void);

/** Append a copy of `polygon`, growing the set if it must.
 *
 * Polygons of fewer than three vertices, and those thinner than
 * `CHOPPER_POLYGON_AREA_TOLERANCE` of their own bounding box, are dropped rather than
 * stored: they carry no phase space and would otherwise accumulate.
 *
 * @return 0 if the set could not grow, 1 otherwise -- including when the polygon was
 *         deliberately dropped, which is not a failure
 */
int chopper_polygon_set_add(chopper_polygon_set * set, const chopper_polygon * polygon);

/** Release what a set allocated and leave an empty set behind. A NULL pointer is fine. */
void chopper_polygon_set_free(chopper_polygon_set * set);

/** Total area, s^2/m -- a plain sum, because the pieces are disjoint. */
double chopper_polygon_set_area(const chopper_polygon_set * set);

/** Replace a set with the part of it that passes one disk.
 *
 * Each polygon contributes one output per turn of each opening it can reach, so the set
 * grows before it shrinks; for a real train most of those turns miss and the survivors
 * collapse quickly. A disk parked open is skipped, having no period to constrain anything,
 * and one parked shut empties the set and says so on stdout, as everywhere else here. A
 * disk open for at least a whole turn is skipped for the same reason as one parked open.
 *
 * `path_spread` is how much further than `chopper.path` a neutron may have travelled to
 * reach this disk -- a guide's path-length spread, in metres. Zero is the straight line
 * every other function in this file assumes. It widens the answer, and is a *support*
 * rather than a distribution: a neutron is passed if some path in `[path, path +
 * path_spread]` would have got it through, with no weighting over which. That is the same
 * bargain `chopper_parameters::aperture` makes for the width of the beam.
 *
 * A spread wide enough to make consecutive turns of one disk overlap would break the
 * disjointness the area and the sampler rely on, and be silently wrong rather than loudly
 * so. This refuses instead: the condition is
 *
 *     path_spread * largest inverse velocity in the set < period - opening
 *
 * and a real guide is orders of magnitude inside it.
 *
 * @param set The set, replaced in place; left empty and owning nothing on failure
 * @param chopper The disk
 * @param path_spread Extra path available, m; 0 for a straight line
 * @return 0 if memory ran out, a polygon overflowed, or the spread would overlap turns;
 *         1 otherwise
 */
int chopper_polygon_set_transmit(chopper_polygon_set * set, chopper_parameters chopper,
                                 double path_spread);

/** The whole train, in beam order, stopping early once nothing is left.
 *
 * The result is the exact transmitted region: `chopper_polygon_set_area` over the area of
 * what was handed in is the fraction of a uniformly drawn ray budget the train passes, and
 * `chopper_polygon_set_inverse_velocity_ranges` is the band list.
 *
 * @param set The source region, replaced by what it transmits
 * @param count How many disks
 * @param choppers The disks
 * @param path_spreads One per disk, or NULL for the straight line to every one
 * @return 0 on failure, as `chopper_polygon_set_transmit`
 */
int chopper_polygon_set_transmit_train(chopper_polygon_set * set, unsigned count,
                                       const chopper_parameters * choppers,
                                       const double * path_spreads);

/** The inverse velocity bands a set covers, sorted and merged.
 *
 * What `chopper_inverse_velocity_windows` is trying to compute. That one projects each
 * disk separately and intersects the projections; this projects the intersection, which is
 * the question that was being asked.
 *
 * @warning The returned `ranges` is allocated here and must be freed at calling scope.
 */
range_set chopper_polygon_set_inverse_velocity_ranges(const chopper_polygon_set * set);

/** Whether a point lies inside a convex polygon, its boundary counting as inside.
 *
 * The test every edge in turn: a point inside a convex polygon is on the same side of all
 * of them. `#pragma acc routine seq`, so a McStas TRACE can call it per ray.
 *
 * @param polygon The polygon; fewer than three vertices contains nothing
 * @param inverse_velocity s/m
 * @param time s, at the source
 * @return 1 if the point is inside or on the boundary, 0 otherwise
 */
int chopper_polygon_contains(const chopper_polygon * polygon,
                             double inverse_velocity, double time);

/** Whether a point lies in any polygon of a set.
 *
 * Linear in the total vertex count, which for a real train is a handful: a rectangle
 * through the six BIFROST disks leaves one polygon of five vertices. Also
 * `#pragma acc routine seq`.
 */
int chopper_polygon_set_contains(const chopper_polygon_set * set,
                                 double inverse_velocity, double time);

/** A direct sampler over a transmitted region.
 *
 * The same job as `chopper_mask_sampler` and the same shape -- three uniform deviates, one
 * binary search, never rejects -- over the exact region rather than a grid approximation
 * of it. Each polygon is fanned into triangles from its first vertex and a point is placed
 * in one of them, so `count` here is triangles, not polygons.
 *
 * `acceptance` is exact: the transmitted area over the area sampled, both known in closed
 * form rather than counted in cells. The grid sampler can only over-estimate it, because a
 * partly covered cell is weighted whole -- on a BIFROST train over a wide source that is
 * 6.7% high on a twelve-million-cell grid, and it improves only as the cell count.
 *
 * Everything the note on `chopper_mask_sampler::acceptance` says about independence still
 * applies: the correction is right only while both coordinates are drawn uniformly and
 * independently of each other and of whatever else the caller samples.
 *
 * @param count Triangles to draw from; 0 means nothing to draw
 * @param acceptance Transmitted area over sampled area, and the factor a ray weight must
 *                   be multiplied by
 * @param cumulative `count` entries increasing to 1, the area-weighted triangle distribution
 * @param origin First vertex of each triangle
 * @param edge_a Second vertex less the first
 * @param edge_b Third vertex less the first
 */
struct chopper_polygon_sampler_struct {
  unsigned count;
  double acceptance;
  double * cumulative;
  chopper_point * origin;
  chopper_point * edge_a;
  chopper_point * edge_b;
};
typedef struct chopper_polygon_sampler_struct chopper_polygon_sampler;

/** Build a sampler over a transmitted region.
 *
 * @param set What the train transmits
 * @param sampled_area The area of the region the caller draws from, s^2/m -- normally
 *                     `chopper_polygon_area` of the source rectangle handed to the train.
 *                     It sets `acceptance` and nothing else.
 * @return A sampler whose arrays are allocated here and must be released with
 *         `chopper_polygon_sampler_free`. An empty set, or a non-positive `sampled_area`,
 *         comes back with `count` 0 and no allocations; drawing from it is a caller error.
 */
chopper_polygon_sampler chopper_polygon_sampler_make(const chopper_polygon_set * set,
                                                     double sampled_area);

/** Zero a sampler so it can be freed or tested before it has been built. */
void chopper_polygon_sampler_empty(chopper_polygon_sampler * sampler);

/** Release what `chopper_polygon_sampler_make` allocated. A NULL pointer is fine. */
void chopper_polygon_sampler_free(chopper_polygon_sampler * sampler);

/** Draw one (inverse velocity, time) pair uniformly from the transmitted region.
 *
 * As with `chopper_mask_sampler_draw`, the deviates are arguments rather than drawn here,
 * so the library needs no generator and a caller inside a McStas TRACE can hand over
 * `rand01()`.
 *
 * @param sampler A sampler with a non-zero `count`
 * @param triangle_deviate A uniform deviate on [0, 1), choosing which triangle
 * @param first_deviate A uniform deviate on [0, 1), placing the point along one edge
 * @param second_deviate A uniform deviate on [0, 1), along the other
 * @param inverse_velocity [out] The drawn inverse velocity, s/m
 * @param time [out] The drawn emission time, s
 *
 * @note The caller still owes the ray weight a factor of `sampler->acceptance`.
 */
void chopper_polygon_sampler_draw(const chopper_polygon_sampler * sampler,
                                  double triangle_deviate, double first_deviate,
                                  double second_deviate,
                                  double * inverse_velocity, double * time);

int chopper_write_mask_to_file(
  const char * directory, const char * filename, const char * extension, const char * path_sep,
  const int * mask, unsigned inverse_velocity_count, unsigned time_count,
  const double * inverse_velocities, const double * times
);

/** Write a transmitted region as JSON.
 *
 * The grid writers above put a picture on a fixed mesh; this writes the region itself, so
 * nothing is quantised and the file is a few hundred bytes rather than a few megabytes.
 * Every number is written with enough digits to read back bit-exact.
 *
 *     {
 *       "chopper_lib_version": "4.2.0",
 *       "inverse_velocity_unit": "s/m",
 *       "time_unit": "s",
 *       "sampled": {"inverse_velocity": [lo, hi], "time": [lo, hi], "area": A},
 *       "transmitted_area": a,
 *       "acceptance": a / A,
 *       "inverse_velocity_bands": [[lo, hi], ...],
 *       "polygons": [{"area": ..., "vertices": [[iv, t], ...]}, ...]
 *     }
 *
 * @param directory Where to write, or NULL
 * @param filename The base name
 * @param extension Appended unless `filename` already ends with it
 * @param path_sep The platform's separator, as a string
 * @param set The transmitted region
 * @param sampled The region the caller drew from, which sets `acceptance`; may be NULL,
 *                and then `sampled` and `acceptance` are written as null
 * @return 1 on success, 0 if the file could not be opened
 */
int chopper_write_polygons_to_file(
  const char * directory, const char * filename, const char * extension, const char * path_sep,
  const chopper_polygon_set * set, const chopper_polygon * sampled
);

int chopper_write_total_to_file(
  const char * directory, const char * filename, const char * extension, const char * path_sep,
  const double * total, unsigned inverse_velocity_count, unsigned time_count,
  const double * inverse_velocities, const double * times
);

#endif //CHOPPER_LIB_CHOPPER_LIB_H
