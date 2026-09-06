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
#define CHOPPER_LIB_VERSION_MINOR 0
#define CHOPPER_LIB_VERSION_PATCH 0
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

/** Sort the limits of a single `range` in place
 *
 * A range is characterized by its minimal and maximal edges. This function ensures they are ordered properly.
 * */
void range_sort(range a);

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
 */
struct chopper_parameters_struct {
  double speed; // rotation frequency in Hz
  double delay; // when the disk point at angle `beam` is on the path, in seconds
  double beam; // from the zero mark to the beam crossing, in degrees
  unsigned edge_count; // number of entries in edges, two per opening
  double * edges; // opening and closing edge of each opening, in degrees from the mark
  double path; // average(?) path length from source to this chopper in meters
};
typedef struct chopper_parameters_struct chopper_parameters;


/** Find the possible inverse velocity window(s) that are admitted by a series of disk choppers
 *
 * @param count The number of disk choppers provided
 * @param choppers The parameters of the disk choppers
 * @param inv_v_min The minimum inverse velocity to be considered -- likely matching a guide cutoff
 * @param inv_v_max The maximum inverse velocity to be considered -- how long before a neutron is no-longer interesting
 * @param latest_emission How long after time-zero can a neutron start its journey, effects minimum inverse velocities
 * @return One or more inverse velocity ranges that can pass through the chopper train as a `range_set`
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

int chopper_write_mask_to_file(
  const char * directory, const char * filename, const char * extension, const char * path_sep,
  const int * mask, unsigned inverse_velocity_count, unsigned time_count,
  const double * inverse_velocities, const double * times
);

int chopper_write_total_to_file(
  const char * directory, const char * filename, const char * extension, const char * path_sep,
  const double * total, unsigned inverse_velocity_count, unsigned time_count,
  const double * inverse_velocities, const double * times
);

#endif //CHOPPER_LIB_CHOPPER_LIB_H
