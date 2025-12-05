#ifndef CHOPPER_LIB_CHOPPER_LIB_H
#define CHOPPER_LIB_CHOPPER_LIB_H
//
// Created by Gregory Tucker, ESS ERIC on 2023-06-01.
//

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


/** The parameters of a single-opening disk chopper.
 *
 * @param speed The rotation speed of the disk, in Hz
 * @param phase The orientation of the disk at 'zero'-time, in degrees
 * @param angle The opening size of the disk, in degrees
 * @param path  The path length from the 'zero'-time source to the disk positon, in meters
 *
 * A multi-opening chopper could be treated as a single-opening chopper if all openings are the same size and
 * are distributed equally around the disk. In such a case the speed parameter of this structure should be the
 * 'opening appearance' frequency, so the rotation speed times the number of equally spaced openings.
 */
struct chopper_parameters_struct {
  double speed; // rotation frequency in Hz
  double phase; // t=0 orientation in degrees
  double angle; // *single* window opening angle in degrees
  double path; // average(?) path length from source to this chopper in meters
};
typedef struct chopper_parameters_struct chopper_parameters;

struct chopper_window_struct {
  double min; // the minimum angle of the window in degrees
  double max; // the maximum angle of the window in degrees
};
typedef struct chopper_window_struct chopper_window;
struct multi_chopper_parameters_struct {
  double speed; // rotation frequency in Hz
  double phase; // t=0 orientation in degrees
  unsigned window_count; // number of windows
  chopper_window * windows; // array of window definitions
  double path; // average(?) path length from source to this chopper in meters
};
typedef struct multi_chopper_parameters_struct multi_chopper_parameters;

multi_chopper_parameters single_to_multi_chopper(chopper_parameters single);

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

/** Create a mask of allowed (inverse_velocity, time) bins based on multi-chopper parameters
 *
 * @param mask [out] An array of integers to be filled with 1 (allowed) or 0 (blocked), of size (inverse_velocity_count-1) * (time_count-1)
 * @param mask_inverse_velocity_count [in] The number of inverse velocity bins in the mask (should be inverse_velocity_count - 1)
 * @param mask_time_count [in] The number of time bins in the mask (should be time_count - 1)
 * @param inverse_velocities [in] An array of inverse velocities (in s/m), of size inverse_velocity_count
 * @param inverse_velocity_count [in] The number of inverse velocities provided
 * @param times [in] An array of times at the source position (in s), of size time_count
 * @param time_count [in] The number of times provided
 * @param choppers [in] An array of multi chopper parameters, of size chopper_count
 * @param chopper_count [in] The number of choppers provided
 * @param grow_mask [in] Expand the allowed regions by this number of bins in each direction
 * @return The number of unmasked (allowed) (inverse_velocity, time) bins
 */
unsigned multi_chopper_inverse_velocity_time_mask(
  int * mask, unsigned mask_inverse_velocity_count, unsigned mask_time_count,
  const double * inverse_velocities, unsigned inverse_velocity_count,
  const double * times, unsigned time_count,
  const multi_chopper_parameters * choppers, unsigned chopper_count,
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
