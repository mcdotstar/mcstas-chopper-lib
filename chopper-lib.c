//
// Created by Gregory Tucker, ESS ERIC on 2023-06-01.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef CHOPPER_LIB_CHOPPER_LIB_H
#include "chopper-lib.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/******************************** range functions ************************************/
void range_sort(range a){
  if (a.maximum < a.minimum){
    const double tmp = a.minimum;
    a.minimum = a.maximum;
    a.maximum = tmp;
  }
}
int classify_range_overlap(const range * a, const range * b){
  //  |A    A|   |B   B| ... or ... |B   B|   |A   A|
  if (a->maximum < b->minimum || b->maximum < a->minimum) return 0;
  // there should be *some* overlap at this point:
  if (a->minimum == b->minimum && a->maximum == b->maximum) return 1; // identical
  // |2| -> B on right end; |3| -> A on right end; (+) -> one inside the other, (-) -> overlapping subregion
  if (a->minimum <= b->minimum && a->maximum >= b->maximum) return 3; // |A|BB|A|
  if (b->minimum <= a->minimum && b->maximum >= a->maximum) return 2; // |B|AA|B|
  if (a->minimum < b->minimum && a->maximum < b->maximum) return -2; // |A|BA|B|
  if (b->minimum < a->minimum && b->maximum < a->maximum) return -3; // |B|AB|A|
  return 0;
}
int compare_ranges(const range * a, const range * b){
  if (a->minimum > b->minimum) return 1;
  if (a->minimum < b->minimum) return -1;
  return 0;
}
// This gateway function is used along with qsort, which handles only void pointers
static int compare_sorted_ranges(const void * ptr_a, const void * ptr_b){
  return compare_ranges((range *) ptr_a, (range *) ptr_b);
}
/******************************** range_set functions ************************************/
range_set range_set_sort(range_set s){
  // An empty set is already sorted, and qsort may not be handed a null pointer even for
  // nothing to sort. A set with no ranges is ordinary here rather than exceptional: it is
  // what a chopper admitting no window intersects with, and what a train that passes
  // nothing carries from there on.
  if (s.count == 0 || s.ranges == NULL) return s;
  // sort all sub-ranges:
  for (unsigned i=0; i<s.count; ++i) range_sort(s.ranges[i]);
  // sort the sub-ranges by minimum
  qsort(s.ranges, s.count, sizeof(*(s.ranges)), compare_sorted_ranges);
  // combine overlapping ranges:
  unsigned overlapping = 0;
  for (unsigned i=1; i<s.count; ++i) if (s.ranges[i-1].maximum >= s.ranges[i].minimum) ++overlapping;
  if (overlapping){
    range_set new_s;
    new_s.count = s.count - overlapping;
    new_s.ranges = calloc(new_s.count, sizeof(range));
    // copy the first element
    new_s.ranges[0].minimum = s.ranges[0].minimum;
    new_s.ranges[0].maximum = s.ranges[0].maximum;

    unsigned copied = 1;
    for (unsigned i=1; i<s.count; ++i) if (s.ranges[i-1].maximum >= s.ranges[i].minimum) {
        // combine the lower bound of the end of the new ranges and the i_th range upper bound over the last new range:
        new_s.ranges[copied-1].maximum = s.ranges[i].maximum;
      } else {
        new_s.ranges[copied].minimum = s.ranges[i].minimum;
        new_s.ranges[copied++].maximum = s.ranges[i].maximum;
      }
    if (copied != s.count - overlapping) printf("Expected to copy %u but copied %u ranges!\n", s.count - overlapping, copied);
    // // free the now-old range_set before we lose its handle ... this is dangerous
    // if (s.ranges) free(s.ranges);
    // recursively re-sort in case we missed overlapping ranges
    return range_set_sort(new_s);
  } else {
    return s;
  }
}

static int range_intersects_ranges(const range r, const range_set rs){
  for (unsigned i=0; i<rs.count; ++i){
    if (classify_range_overlap(&r, rs.ranges + i) != 0) return 1;
  }
  return 0;
}

range_set range_intersection(const range_set ain, const range_set bin){
  range_set a = range_set_sort(ain);
  range_set b = range_set_sort(bin);
  range_set out;
  out.count = 0;
  unsigned i=0, j=0;
  while (i < a.count && j < b.count){
    switch (abs(classify_range_overlap(a.ranges + i, b.ranges + j))) {
      case 3: out.count++; ++j; break; // -3: |B|AB|A|, 3: |A|BB|A|, both increment B
      case 2: out.count++; ++i; break; // -2: |A|BA|B|, 2: |B|AA|B|, both increment A
      case 1: out.count++; ++i; ++j; break; // identical, increment both
      default: {
        // no overlap, so increment the range with the smaller minimum edge
        int comp = compare_ranges(a.ranges + i, b.ranges + j);
        if (comp < 0) ++i;
        if (comp > 0) ++j;
        if (comp == 0) {
          printf("This should not be possible");
          exit(-1);
        }
      }
    }
  }
  // we now know *how many* output sub-ranges there will be, and can allocate the output structure
  out.ranges = calloc(out.count, sizeof(range));
  // and can go through again actually assigning outputs:
  i = 0;
  j = 0;
  unsigned k=0;
  while (i < a.count && j < b.count){
    switch (classify_range_overlap(a.ranges + i, b.ranges + j)) {
      case -3: { // |B|AB|A|; keep (A[i]min, B[j]max) increment j since A extends to higher value
        out.ranges[k].minimum = a.ranges[i].minimum;
        out.ranges[k].maximum = b.ranges[j].maximum;
        ++j; ++k;
        break;
      }
      case -2: { //|A|BA|B|; keep (B[j]min, A[i]max) increment i since B extends to higher value
        out.ranges[k].minimum = b.ranges[j].minimum;
        out.ranges[k].maximum = a.ranges[i].maximum;
        ++i; ++k;
        break;
      }
      case 1: { // identical, keep either a[i] or b[j] and increment both i & j
        out.ranges[k].minimum = a.ranges[i].minimum;
        out.ranges[k].maximum = a.ranges[i].maximum;
        ++i; ++j; ++k;
        break;
      }
      case 2: { // |B|AA|B|; keep a[i] and increment i
        out.ranges[k].minimum = a.ranges[i].minimum;
        out.ranges[k].maximum = a.ranges[i].maximum;
        ++i; ++k;
        break;
      }
      case 3: { // |A|BB|A|; keep b[j] and increment j
        out.ranges[k].minimum = b.ranges[j].minimum;
        out.ranges[k].maximum = b.ranges[j].maximum;
        ++j; ++k;
        break;
      }
      default: {
        if (compare_ranges(a.ranges + i, b.ranges + j) < 0) ++i; else ++j;
      }
    }
  }
  // memory management: (did the sort function make a new range_set?)
  if (a.ranges != NULL && a.ranges != ain.ranges) free(a.ranges);
  if (b.ranges != NULL && b.ranges != bin.ranges) free(b.ranges);
  return out;
}

/****************** chopper train functionality ************************/
/** When a disk edge at angle `a` is on the beam, ignoring which rotation.
 *
 * The angular term is negated because `edges` increase the way NXdisk_chopper measures
 * them, so a turning disk carries a larger angle *towards* the beam rather than away
 * from it; and it keeps the sign of `speed`, so reversing the disk places the opening on
 * the other side of `delay`. `beam` is where the beam crosses the disk, measured from
 * the same mark the edges are.
 */
static double chopper_edge_time(const chopper_parameters chopper, const double a) {
  return chopper.delay + (chopper.beam - a) / 360.0 / chopper.speed;
}

/** How long before and after those times a beam of finite width is still in the opening.
 *
 * An opening is angular and a beam is not. A neutron crossing the disk to one side of the
 * beam centre meets an edge before or after one crossing at the centre does, by the angle
 * between them over the rate the disk turns. Half the aperture either side of every edge
 * time is the whole of it, and it keeps the sign of nothing: a disk turning backwards is
 * as wide as one turning forwards, so the period, not the signed speed, sets the scale.
 */
static double chopper_aperture_time(const chopper_parameters chopper) {
  if (chopper.aperture <= 0 || chopper.speed == 0) return 0.0;
  return chopper.aperture / 2.0 / 360.0 / fabs(chopper.speed);
}

/** Whether a disk that is not turning stands open on the beam.
 *
 * A parked disk is open or shut for good. The beam crosses it at `beam`, the openings
 * span `edges` in the same frame, and with no speed there is no period and no delay to
 * apply: either that angle is inside an opening or it is on the solid part of the disk.
 * Which one decides whether the disk drops out of a calculation or empties it.
 */
int chopper_parked_is_open(const chopper_parameters chopper) {
  if (chopper.edge_count < 2 || chopper.edges == NULL) return 0;
  for (unsigned i = 0; i + 1 < chopper.edge_count; i += 2) {
    const double width = chopper.edges[i + 1] - chopper.edges[i];
    // measured from this opening's first edge, so an opening written across the mark --
    // {350, 370} -- needs no special case, and neither does a negative `beam`
    double from_edge = fmod(chopper.beam - chopper.edges[i], 360.0);
    if (from_edge < 0) from_edge += 360.0;
    if (from_edge < width) return 1;
  }
  return 0;
}

/** Name the disk that emptied the answer.
 *
 * A disk parked shut is a beam stop, and the answer either function can give for one is
 * empty -- the same answer an over-constrained train gives, and the same one a disk with
 * no openings gives. An empty range set says nothing about which disk emptied it, and a
 * disk parked shut is nearly always one left mis-set rather than the question being
 * asked, so say which one it was and why.
 */
static void chopper_report_parked_shut(const chopper_parameters chopper, const unsigned index) {
  printf("chopper-lib: nothing gets through chopper %u, parked with the beam at %g "
         "degrees, where the disk is solid; the train admits nothing.\n",
         index, chopper.beam);
}

range_set chopper_inverse_velocity_windows(const unsigned count, const chopper_parameters * choppers,
                                           const double inv_v_min, const double inv_v_max,
                                           const double latest_emission){
  range_set limits;
  limits.count = 1;
  limits.ranges = calloc(1, sizeof(range));
  limits.ranges[0].minimum = inv_v_min;
  limits.ranges[0].maximum = inv_v_max;

  for (unsigned i=0; i<count && limits.count; ++i) {
    // A disk that is not turning has no period, so it admits every time or none. Parked
    // open it constrains nothing and drops out; parked shut it is a beam stop, exactly
    // like a disk with no openings at all, and no other disk can undo that.
    if (choppers[i].speed == 0.0) {
      if (chopper_parked_is_open(choppers[i])) continue;
      chopper_report_parked_shut(choppers[i], i);
      if (limits.ranges) free(limits.ranges);
      limits.count = 0;
      limits.ranges = NULL;
      break;
    }
    // the period of the chopper is a positive time
    const double tau = 1.0 / fabs(choppers[i].speed);
    const double path = choppers[i].path;
    const double aperture = chopper_aperture_time(choppers[i]);
    const unsigned opening_count = choppers[i].edge_count / 2;
    // allocate open and close time arrays for the openings to avoid the same calculation twice
    double * t_open = (double *) calloc(opening_count, sizeof(double));
    double * t_close = (double *) calloc(opening_count, sizeof(double));
    // find the overall minimum and maximum total rotations that cover the inverse velocity range
    int first=1, n_min=0, n_max=0;
    for (unsigned opening=0; opening<opening_count; ++opening) {
      // Which edge of an opening reaches the beam first depends on which way the disk
      // turns, so place both and order the pair afterwards rather than by name.
      t_open[opening] = chopper_edge_time(choppers[i], choppers[i].edges[2 * opening]);
      t_close[opening] = chopper_edge_time(choppers[i], choppers[i].edges[2 * opening + 1]);
      if (t_close[opening] < t_open[opening]) {
        const double tmp = t_open[opening];
        t_open[opening] = t_close[opening];
        t_close[opening] = tmp;
      }
      // and a beam of finite width reaches the opening early and leaves it late
      t_open[opening] -= aperture;
      t_close[opening] += aperture;
      const int n_j_min = (int) floor((path * inv_v_min - t_open[opening]) / tau);
      const int n_j_max = (int) ceil((path * inv_v_max - t_close[opening]) / tau);
      if (first || n_j_min < n_min) n_min = n_j_min;
      if (first || n_j_max > n_max) n_max = n_j_max;
      if (first) first = 0;
    }
    // collect the ranges for each of the n in (n_min, n_max) into a set:
    const unsigned rotation_count = (unsigned)(n_max - n_min + 1);
    range_set ith;
    // every rotation can contribute one range per opening, so there has to be room for all of them
    ith.count = rotation_count * opening_count;
    ith.ranges = ith.count ? calloc(ith.count, sizeof(range)) : NULL;
    if (ith.count && ith.ranges == NULL) {
      printf("Out of memory\n");
      exit(-1);
    }
    unsigned c=0;
    for (unsigned j=0; j < rotation_count; ++j) {
      const double n_tau = tau * (double) (n_min + (int) j);
      for (unsigned opening=0; opening<opening_count; ++opening) {
        // let the minimum 1/v come from the *end* of the pulse:
        double wiv_min = (t_open[opening] + n_tau - latest_emission) / path;
        double wiv_max = (t_close[opening] + n_tau) / path;
        // clamp to provided limits
        wiv_min = wiv_min < inv_v_min ? inv_v_min : wiv_min > inv_v_max ? inv_v_max : wiv_min;
        wiv_max = wiv_max < inv_v_min ? inv_v_min : wiv_max > inv_v_max ? inv_v_max : wiv_max;
        // and insert if this is a valid range:
        if (wiv_min < wiv_max) {
          ith.ranges[c].minimum = wiv_min;
          ith.ranges[c++].maximum = wiv_max;
        }
      }
    }
    ith.count = c;
    // clean up allocated opening and closing times
    free(t_open);
    free(t_close);
    // find the intersection of this chopper with the running set
    const range_set new_limits = range_intersection(limits, ith);
    // clean-up allocated memory, ensuring limits or ith is not erased if transferred to new_limits:
    if ((!new_limits.count || new_limits.ranges != limits.ranges) && limits.ranges)  free(limits.ranges);
    if ((!new_limits.count || new_limits.ranges != ith.ranges) && ith.ranges) free(ith.ranges);
    limits = new_limits;
  }
  return limits;
}

unsigned chopper_inverse_velocity_limits(double * lower, double * upper,
                                         const unsigned count, const chopper_parameters * choppers,
                                         const double inv_v_min, const double inv_v_max,
                                         const double latest_emission){
  const range_set limits = chopper_inverse_velocity_windows(count, choppers, inv_v_min, inv_v_max, latest_emission);
  if (limits.count){
    *lower = limits.ranges[0].minimum;
    *upper = limits.ranges[limits.count-1].maximum;
  }
  if (limits.ranges) free(limits.ranges);
  return limits.count;
}

// V2K, K2V and PI come from whatever is compiling this. McStas defines all three in
// its runtime, and this file is copied verbatim into every instrument that
// %includes it, so defining them here would put a second definition of PI in the
// generated C. Every other build passes them in instead -- CMake does it with
// CHOPPER_LIB_DEFINITIONS, which the README spells out for a consumer that compiles
// this source into a target of its own rather than linking the library.
#if !defined(V2K) || !defined(K2V) || !defined(PI)
#error "chopper-lib.c needs V2K, K2V and PI defined; build it with CHOPPER_LIB_DEFINITIONS or inside McStas"
#endif

unsigned chopper_wavelength_limits(double * lower, double * upper,
                                   const unsigned count, const chopper_parameters * choppers,
                                   const double lambda_min, const double lambda_max, const double latest_emission){
  const unsigned windows = chopper_inverse_velocity_limits(
    lower, upper, count, choppers, lambda_min * V2K / 2 / PI, lambda_max * V2K / 2 / PI, latest_emission);
  if (windows) {
    *lower *= K2V * 2 * PI;
    *upper *= K2V * 2 * PI;
  }
  return windows;
}

static int_range chopper_rotation_limits(const chopper_parameters chopper, const range time_range) {
  int_range rotations = {.minimum = 1, .maximum = -1};
  if (chopper.edge_count < 2 || chopper.edges == NULL) {
    return rotations;
  }
  // The edges increase, so the disk's angular extent is its first and its last -- which
  // is the one simplification the flat array buys outright over a list of pairs.
  const range edges_range = {.minimum = chopper.edges[0],
                             .maximum = chopper.edges[chopper.edge_count - 1]};
  // An edge at angle a is on the beam at chopper_edge_time(a), and every period
  // thereafter, so the rotation index of a time t is (t - edge_time) / tau. Both terms
  // must be in rotations: dividing the elapsed time by the period is what converts it,
  // and the period is positive however the disk turns.
  const double tau = 1.0 / fabs(chopper.speed);
  const double t0 = chopper.delay;
  // Convert the disk's angular extent to fractional rotations, keeping the pair ordered:
  // reversing the disk exchanges which edge is the earlier one.
  double lowest = (chopper_edge_time(chopper, edges_range.maximum) - t0) / tau;
  double highest = (chopper_edge_time(chopper, edges_range.minimum) - t0) / tau;
  if (lowest > highest) {
    const double tmp = lowest;
    lowest = highest;
    highest = tmp;
  }
  // a beam of finite width opens the disk earlier and closes it later, in the same
  // rotations, so the extent it can cover is that much wider at both ends
  const double aperture = chopper_aperture_time(chopper) / tau;
  lowest -= aperture;
  highest += aperture;
  // find the number of rotations needed to place the latest angle _before_ the earliest time:
  rotations.minimum = (int) floor((time_range.minimum - t0) / tau - highest) - 1;
  // and the number of rotations needed to place the earliest angle _after_ the latest time:
  rotations.maximum = (int) ceil((time_range.maximum - t0) / tau - lowest) + 1;
  return rotations;
}

unsigned chopper_inverse_velocity_time_mask(
  int *mask, const unsigned mask_inverse_velocity_count, const unsigned mask_time_count,
  const double *inverse_velocities, const unsigned inverse_velocity_count,
  const double *times, const unsigned time_count,
  const chopper_parameters *choppers, const unsigned chopper_count,
  const int grow_mask
  ) {
  unsigned allowed_bins = 0;

  if (mask_inverse_velocity_count != inverse_velocity_count - 1 || mask_time_count != time_count - 1) {
    printf("Mask dimensions do not match the provided inverse velocity and time arrays\n");
    exit(-1);
  }

  // Initialize the edges for bin edge checks
  const unsigned inverse_velocity_edges_count = inverse_velocity_count * mask_time_count;
  const unsigned time_edges_count = mask_inverse_velocity_count * time_count;
  int * inverse_velocity_edges = calloc(inverse_velocity_edges_count, sizeof(int));
  if (inverse_velocity_edges == NULL) {
    printf("Out of memory\n");
    exit(-1);
  }
  int * time_edges = calloc(time_edges_count, sizeof(int));
  if (time_edges == NULL) {
    printf("Out of memory\n");
    exit(-1);
  }
  for (unsigned i = 0; i < inverse_velocity_edges_count; ++i) inverse_velocity_edges[i] = 1;
  for (unsigned i = 0; i < time_edges_count; ++i) time_edges[i] = 1;

  for (unsigned ci = 0; ci < chopper_count; ++ci) {
    // A stationary disk has no period to speak of; the window functions skip one parked
    // open rather than dividing by zero, so do the same here. One parked shut takes the
    // beam stop path below, which is what it is.
    if (choppers[ci].speed == 0.0) {
      if (chopper_parked_is_open(choppers[ci])) continue;
      chopper_report_parked_shut(choppers[ci], ci);
      memset(inverse_velocity_edges, 0, inverse_velocity_edges_count * sizeof(int));
      memset(time_edges, 0, time_edges_count * sizeof(int));
      break; // nothing any other chopper does can let a neutron back through
    }
    // A disk with no openings is a beam stop rather than an absent chopper, which is what
    // the window functions report for the same chopper: their range set comes back empty.
    if (choppers[ci].edge_count < 2 || choppers[ci].edges == NULL) {
      memset(inverse_velocity_edges, 0, inverse_velocity_edges_count * sizeof(int));
      memset(time_edges, 0, time_edges_count * sizeof(int));
      break; // nothing any other chopper does can let a neutron back through
    }
    const double tau = 1.0 / fabs(choppers[ci].speed);
    const double aperture = chopper_aperture_time(choppers[ci]);
    const range time_range = {
      .minimum = times[0] + choppers[ci].path * inverse_velocities[0],
      .maximum = times[time_count - 1] + choppers[ci].path * inverse_velocities[inverse_velocity_count - 1]
    };
    const int_range rotations = chopper_rotation_limits(choppers[ci], time_range);
    if (rotations.maximum < rotations.minimum) continue; // No possible rotations
    // build the ranges of allowed times for this chopper
    range_set allowed_times;
    const unsigned opening_count = choppers[ci].edge_count / 2;
    allowed_times.count = (rotations.maximum - rotations.minimum + 1) * opening_count;
    allowed_times.ranges = calloc(allowed_times.count, sizeof(range));
    unsigned c = 0;
    for (int n = rotations.minimum; n <= rotations.maximum; ++n) {
      for (unsigned w = 0; w < opening_count; ++w) {
        // Which edge of an opening reaches the beam first depends on which way the disk
        // turns, so place both and order the pair afterwards.
        const double a = (double) n * tau + chopper_edge_time(choppers[ci], choppers[ci].edges[2 * w]);
        const double b = (double) n * tau + chopper_edge_time(choppers[ci], choppers[ci].edges[2 * w + 1]);
        // widened by the beam's own width, which is why this is not the same as growing
        // the finished mask: it opens the window in time and leaves the velocities alone
        allowed_times.ranges[c].minimum = (a < b ? a : b) - aperture;
        allowed_times.ranges[c++].maximum = (a < b ? b : a) + aperture;
      }
    }
    // Now check each bin edge against the allowed times
    for (unsigned ti = 0; ti < time_count; ++ti) {
      for (unsigned vi = 0; vi < mask_inverse_velocity_count; ++ vi) {
        if (time_edges[ti * mask_inverse_velocity_count + vi] == 1) {
          const range edge= {
            .minimum = times[ti] + choppers[ci].path * inverse_velocities[vi],
            .maximum = times[ti] + choppers[ci].path * inverse_velocities[vi + 1]
          };
          if (range_intersects_ranges(edge, allowed_times) == 0) {
            time_edges[ti * mask_inverse_velocity_count + vi] = 0;
          }
        }
      }
    }
    for (unsigned vi = 0; vi < inverse_velocity_count; ++vi) {
      for (unsigned ti = 0; ti < mask_time_count; ++ti) {
        if (inverse_velocity_edges[vi * mask_time_count + ti] == 1) {
          const range edge = {
            .minimum = times[ti] + choppers[ci].path * inverse_velocities[vi],
            .maximum = times[ti + 1] + choppers[ci].path * inverse_velocities[vi]
          };
          if (range_intersects_ranges(edge, allowed_times) == 0) {
            inverse_velocity_edges[vi * mask_time_count + ti] = 0;
          }
        }
      }
    }
    // end of this chopper
    free(allowed_times.ranges);
  }
  // go through the edges and mark bins which are disallowed
  for (unsigned ti = 0; ti < mask_time_count; ++ti) {
    for (unsigned vi = 0; vi < mask_inverse_velocity_count; ++vi) {
      const unsigned index = ti * mask_inverse_velocity_count + vi;
      if (inverse_velocity_edges[vi * mask_time_count + ti] == 0 && inverse_velocity_edges[(vi + 1) * mask_time_count + ti] == 0 &&
          time_edges[ti * mask_inverse_velocity_count + vi] == 0 && time_edges[(ti + 1) * mask_inverse_velocity_count + vi] == 0) {
        if (mask != NULL) mask[index] = CHOPPER_MASK_EXCLUDED;
      } else {
        if (mask != NULL) mask[index] = CHOPPER_MASK_INCLUDED;
        ++allowed_bins;
      }
    }
  }

  // Optionally grow the mask in each direction by the number of grow_mask specified
  if (grow_mask && mask != NULL) {
    for (unsigned ti = 0; ti < mask_time_count; ++ti) {
      for (unsigned vi = 0; vi < mask_inverse_velocity_count; ++vi) {
        const unsigned index = ti * mask_inverse_velocity_count + vi;
        if (mask[index] == CHOPPER_MASK_EXCLUDED) {
          int surround = 0;
          // Check the surrounding bins
          for (int dt = -grow_mask; dt <= grow_mask; ++dt) {
            for (int dv = -grow_mask; dv <= grow_mask; ++dv) {
              if (dt == 0 && dv == 0) continue; // Skip the center
              const int nti = (int)ti + dt;
              const int nvi = (int)vi + dv;
              if (nti >= 0 && nti < (int)mask_time_count && nvi >= 0 && nvi < (int)mask_inverse_velocity_count) {
                if (mask[nti * mask_inverse_velocity_count + nvi] == CHOPPER_MASK_INCLUDED) {
                  surround = 1;
                  break;
                }
              }
            }
          }
          if (surround) {
            mask[index] = CHOPPER_MASK_GROWN;
            ++allowed_bins;
          }
        }
      }
    }
    // Now ensure all grown bins are set to allowed (1)
    for (unsigned ti = 0; ti < mask_time_count; ++ti) {
      for (unsigned vi = 0; vi < mask_inverse_velocity_count; ++vi) {
        const unsigned index = ti * mask_inverse_velocity_count + vi;
        if (mask[index] == CHOPPER_MASK_GROWN) {
          mask[index] = CHOPPER_MASK_INCLUDED;
        }
      }
    }
  }

  free(inverse_velocity_edges);
  free(time_edges);
  return allowed_bins;
}


double chopper_unmasked_probability(
  const double * signal, const int * mask, const unsigned mask_inverse_velocity_count, const unsigned mask_time_count
) {
  double unmasked_signal = 0.0;
  double total_signal = 0.0;
  for (unsigned ti = 0; ti < mask_time_count; ++ti) {
    for (unsigned vi = 0; vi < mask_inverse_velocity_count; ++vi) {
      const unsigned index = ti * mask_inverse_velocity_count + vi;
      const double bin_signal = signal[index];
      total_signal += bin_signal;
      if (mask[index] == CHOPPER_MASK_INCLUDED) {
        unmasked_signal += bin_signal;
      }
    }
  }
  return total_signal ? unmasked_signal / total_signal : 0.0;
}

/*************************** mask sampling ******************************************/

/* How much of [low, high] lies inside [limit_low, limit_high]; never negative. */
static double chopper_clipped_width(
  const double low, const double high, const double limit_low, const double limit_high
) {
  const double lo = low > limit_low ? low : limit_low;
  const double hi = high < limit_high ? high : limit_high;
  return hi > lo ? hi - lo : 0.0;
}

void chopper_mask_sampler_empty(chopper_mask_sampler * sampler) {
  if (sampler == NULL) return;
  sampler->count = 0;
  sampler->acceptance = 0.0;
  sampler->cumulative = NULL;
  sampler->inverse_velocity_low = NULL;
  sampler->inverse_velocity_width = NULL;
  sampler->time_low = NULL;
  sampler->time_width = NULL;
}

chopper_mask_sampler chopper_mask_sampler_make(
  const int * mask, const unsigned mask_inverse_velocity_count, const unsigned mask_time_count,
  const double * inverse_velocities, const double * times,
  const double inverse_velocity_minimum, const double inverse_velocity_range,
  const double time_minimum, const double time_range
) {
  chopper_mask_sampler sampler;
  chopper_mask_sampler_empty(&sampler);

  if (mask == NULL || inverse_velocities == NULL || times == NULL) return sampler;
  if (mask_inverse_velocity_count == 0 || mask_time_count == 0) return sampler;
  if (inverse_velocity_range <= 0.0 || time_range <= 0.0) {
    printf("A mask sampler needs a region to sample: given %g s/m by %g s\n",
           inverse_velocity_range, time_range);
    return sampler;
  }

  const double inverse_velocity_maximum = inverse_velocity_minimum + inverse_velocity_range;
  const double time_maximum = time_minimum + time_range;

  /* A cell partly outside the sampled region contributes only the part inside it, and one
   * wholly outside contributes nothing and is left out of the sampler entirely. */
  unsigned allowed = 0;
  for (unsigned ti = 0; ti < mask_time_count; ++ti) {
    const double dt = chopper_clipped_width(times[ti], times[ti + 1], time_minimum, time_maximum);
    if (dt <= 0.0) continue;
    for (unsigned vi = 0; vi < mask_inverse_velocity_count; ++vi) {
      if (mask[ti * mask_inverse_velocity_count + vi] == CHOPPER_MASK_EXCLUDED) continue;
      if (chopper_clipped_width(inverse_velocities[vi], inverse_velocities[vi + 1],
                                inverse_velocity_minimum, inverse_velocity_maximum) > 0.0) {
        ++allowed;
      }
    }
  }
  if (allowed == 0) return sampler;

  sampler.cumulative = calloc(allowed, sizeof(double));
  sampler.inverse_velocity_low = calloc(allowed, sizeof(double));
  sampler.inverse_velocity_width = calloc(allowed, sizeof(double));
  sampler.time_low = calloc(allowed, sizeof(double));
  sampler.time_width = calloc(allowed, sizeof(double));
  if (sampler.cumulative == NULL || sampler.inverse_velocity_low == NULL
      || sampler.inverse_velocity_width == NULL || sampler.time_low == NULL
      || sampler.time_width == NULL) {
    printf("Out of memory building a mask sampler over %u cells\n", allowed);
    chopper_mask_sampler_free(&sampler);
    return sampler;
  }

  double area = 0.0;
  unsigned c = 0;
  for (unsigned ti = 0; ti < mask_time_count; ++ti) {
    const double t_low = times[ti] > time_minimum ? times[ti] : time_minimum;
    const double dt = chopper_clipped_width(times[ti], times[ti + 1], time_minimum, time_maximum);
    if (dt <= 0.0) continue;
    for (unsigned vi = 0; vi < mask_inverse_velocity_count; ++vi) {
      if (mask[ti * mask_inverse_velocity_count + vi] == CHOPPER_MASK_EXCLUDED) continue;
      const double dv = chopper_clipped_width(inverse_velocities[vi], inverse_velocities[vi + 1],
                                              inverse_velocity_minimum, inverse_velocity_maximum);
      if (dv <= 0.0) continue;
      sampler.inverse_velocity_low[c] = inverse_velocities[vi] > inverse_velocity_minimum
                                      ? inverse_velocities[vi] : inverse_velocity_minimum;
      sampler.inverse_velocity_width[c] = dv;
      sampler.time_low[c] = t_low;
      sampler.time_width[c] = dt;
      area += dv * dt;
      sampler.cumulative[c] = area;
      ++c;
    }
  }
  sampler.count = c;
  sampler.acceptance = area / (inverse_velocity_range * time_range);
  /* Normalise the running area into a distribution, and pin the last entry rather than
   * leave a draw above it to fall off the end of a binary search. */
  for (unsigned i = 0; i < sampler.count; ++i) sampler.cumulative[i] /= area;
  sampler.cumulative[sampler.count - 1] = 1.0;

  return sampler;
}

void chopper_mask_sampler_free(chopper_mask_sampler * sampler) {
  if (sampler == NULL) return;
  if (sampler->cumulative) free(sampler->cumulative);
  if (sampler->inverse_velocity_low) free(sampler->inverse_velocity_low);
  if (sampler->inverse_velocity_width) free(sampler->inverse_velocity_width);
  if (sampler->time_low) free(sampler->time_low);
  if (sampler->time_width) free(sampler->time_width);
  chopper_mask_sampler_empty(sampler);
}

/* Unconditional, as McCode's own libraries write it. Whether this is callable from a kernel
 * must not rest on the compiler happening to define __GNUC__, which NVHPC does only for
 * glibc's benefit. MSVC warns on the unknown pragma; the generated instrument silences 4068
 * the same way. */
#ifdef _MSC_EXTENSIONS
#pragma warning(disable: 4068)
#endif
#pragma acc routine seq
void chopper_mask_sampler_draw(
  const chopper_mask_sampler * sampler,
  const double cell_deviate, const double inverse_velocity_deviate, const double time_deviate,
  double * inverse_velocity, double * time
) {
  if (sampler == NULL || sampler->count == 0) return;
  /* The first cell whose cumulative share reaches the deviate. */
  unsigned low = 0, high = sampler->count - 1;
  while (low < high) {
    const unsigned mid = low + (high - low) / 2;
    if (sampler->cumulative[mid] <= cell_deviate) low = mid + 1; else high = mid;
  }
  if (inverse_velocity != NULL) {
    *inverse_velocity = sampler->inverse_velocity_low[low]
                      + sampler->inverse_velocity_width[low] * inverse_velocity_deviate;
  }
  if (time != NULL) {
    *time = sampler->time_low[low] + sampler->time_width[low] * time_deviate;
  }
}

/*************************** transmitted phase space *********************************/

chopper_polygon chopper_polygon_rectangle(const double inverse_velocity_minimum,
                                          const double inverse_velocity_range,
                                          const double time_minimum, const double time_range) {
  chopper_polygon polygon;
  polygon.count = 0;
  if (inverse_velocity_range <= 0.0 || time_range <= 0.0) return polygon;
  const double a0 = inverse_velocity_minimum, a1 = inverse_velocity_minimum + inverse_velocity_range;
  const double t0 = time_minimum, t1 = time_minimum + time_range;
  polygon.vertex[0].inverse_velocity = a0; polygon.vertex[0].time = t0;
  polygon.vertex[1].inverse_velocity = a1; polygon.vertex[1].time = t0;
  polygon.vertex[2].inverse_velocity = a1; polygon.vertex[2].time = t1;
  polygon.vertex[3].inverse_velocity = a0; polygon.vertex[3].time = t1;
  polygon.count = 4;
  return polygon;
}

double chopper_polygon_area(const chopper_polygon * polygon) {
  if (polygon == NULL || polygon->count < 3) return 0.0;
  double twice = 0.0;
  for (unsigned i = 0; i < polygon->count; ++i) {
    const chopper_point p = polygon->vertex[i];
    const chopper_point q = polygon->vertex[(i + 1) % polygon->count];
    twice += p.inverse_velocity * q.time - q.inverse_velocity * p.time;
  }
  return fabs(twice) / 2.0;
}

void chopper_polygon_extent(const chopper_polygon * polygon, const double alpha,
                            const double beta, double * lower, double * upper) {
  if (polygon == NULL || polygon->count == 0) return;
  double lo = alpha * polygon->vertex[0].inverse_velocity + beta * polygon->vertex[0].time;
  double hi = lo;
  for (unsigned i = 1; i < polygon->count; ++i) {
    const double v = alpha * polygon->vertex[i].inverse_velocity + beta * polygon->vertex[i].time;
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  }
  if (lower != NULL) *lower = lo;
  if (upper != NULL) *upper = hi;
}

int chopper_polygon_clip_halfplane(chopper_polygon * polygon, const double alpha,
                                   const double beta, const double c) {
  if (polygon == NULL) return 1;
  if (polygon->count < 3) { polygon->count = 0; return 1; }
  chopper_polygon out;
  out.count = 0;
  for (unsigned i = 0; i < polygon->count; ++i) {
    const chopper_point p = polygon->vertex[i];
    const chopper_point q = polygon->vertex[(i + 1) % polygon->count];
    const double dp = alpha * p.inverse_velocity + beta * p.time - c;
    const double dq = alpha * q.inverse_velocity + beta * q.time - c;
    if (dp <= 0.0) {
      if (out.count >= CHOPPER_POLYGON_MAX_VERTICES) return 0;
      out.vertex[out.count++] = p;
    }
    /* The crossing point is a vertex of the result only when the two ends are strictly on
     * opposite sides. A vertex lying *in* the line is already emitted by the test above --
     * counting it as a crossing too would add it twice, and duplicates eat the vertex
     * budget that CHOPPER_POLYGON_MAX_VERTICES is sized against. */
    if ((dp < 0.0 && dq > 0.0) || (dp > 0.0 && dq < 0.0)) {
      if (out.count >= CHOPPER_POLYGON_MAX_VERTICES) return 0;
      const double f = dp / (dp - dq);
      out.vertex[out.count].inverse_velocity =
          p.inverse_velocity + f * (q.inverse_velocity - p.inverse_velocity);
      out.vertex[out.count].time = p.time + f * (q.time - p.time);
      ++out.count;
    }
  }
  if (out.count < 3) out.count = 0;
  *polygon = out;
  return 1;
}

int chopper_polygon_clip_wedge(chopper_polygon * polygon, const double shortest_path,
                               const double longest_path, const double beta,
                               const double lower, const double upper) {
  /* The earliest a neutron can arrive must not be after the window shuts, and the latest
   * must not be before it opens. */
  if (!chopper_polygon_clip_halfplane(polygon, shortest_path, beta, upper)) return 0;
  if (polygon->count == 0) return 1;
  if (!chopper_polygon_clip_halfplane(polygon, -longest_path, -beta, -lower)) return 0;
  return 1;
}

chopper_polygon_set chopper_polygon_set_empty(void) {
  chopper_polygon_set set;
  set.count = 0;
  set.capacity = 0;
  set.polygon = NULL;
  return set;
}

/** Whether a polygon carries no phase space worth keeping.
 *
 * Against its own bounding box, so the test means the same thing whatever units the
 * caller works in. See CHOPPER_POLYGON_AREA_TOLERANCE.
 */
static int chopper_polygon_is_negligible(const chopper_polygon * polygon) {
  if (polygon == NULL || polygon->count < 3) return 1;
  double a_lo = 0, a_hi = 0, t_lo = 0, t_hi = 0;
  chopper_polygon_extent(polygon, 1.0, 0.0, &a_lo, &a_hi);
  chopper_polygon_extent(polygon, 0.0, 1.0, &t_lo, &t_hi);
  const double box = (a_hi - a_lo) * (t_hi - t_lo);
  if (!(box > 0.0)) return 1;
  return chopper_polygon_area(polygon) <= CHOPPER_POLYGON_AREA_TOLERANCE * box;
}

int chopper_polygon_set_add(chopper_polygon_set * set, const chopper_polygon * polygon) {
  if (set == NULL) return 0;
  if (chopper_polygon_is_negligible(polygon)) return 1;
  if (set->count == set->capacity) {
    const unsigned capacity = set->capacity ? 2 * set->capacity : 8;
    chopper_polygon * grown = (chopper_polygon *) realloc(set->polygon,
                                                          capacity * sizeof(chopper_polygon));
    if (grown == NULL) {
      printf("Out of memory holding %u transmitted polygons\n", capacity);
      return 0;
    }
    set->polygon = grown;
    set->capacity = capacity;
  }
  set->polygon[set->count++] = *polygon;
  return 1;
}

void chopper_polygon_set_free(chopper_polygon_set * set) {
  if (set == NULL) return;
  if (set->polygon) free(set->polygon);
  *set = chopper_polygon_set_empty();
}

double chopper_polygon_set_area(const chopper_polygon_set * set) {
  if (set == NULL) return 0.0;
  double total = 0.0;
  for (unsigned i = 0; i < set->count; ++i) total += chopper_polygon_area(&set->polygon[i]);
  return total;
}

/** The two edge times of one opening, ordered, widened by the beam's own width. */
static void chopper_opening_times(const chopper_parameters chopper, const unsigned opening,
                                  double * lower, double * upper) {
  double lo = chopper_edge_time(chopper, chopper.edges[2 * opening]);
  double hi = chopper_edge_time(chopper, chopper.edges[2 * opening + 1]);
  if (hi < lo) { const double swap = lo; lo = hi; hi = swap; }
  const double aperture = chopper_aperture_time(chopper);
  *lower = lo - aperture;
  *upper = hi + aperture;
}

int chopper_polygon_set_transmit(chopper_polygon_set * set, const chopper_parameters chopper,
                                 const double path_spread) {
  if (set == NULL) return 0;
  if (chopper.edge_count < 2 || chopper.edges == NULL) {
    /* No openings at all is a beam stop, the same as one parked shut. */
    chopper_polygon_set_free(set);
    return 1;
  }
  if (chopper.speed == 0.0) {
    if (chopper_parked_is_open(chopper)) return 1;
    printf("chopper-lib: nothing gets through a chopper parked with the beam at %g "
           "degrees, where the disk is solid; the train transmits nothing.\n", chopper.beam);
    chopper_polygon_set_free(set);
    return 1;
  }
  if (path_spread < 0.0) {
    printf("A path spread is a length, so it has no sign; given %g m\n", path_spread);
    chopper_polygon_set_free(set);
    return 0;
  }

  const double tau = 1.0 / fabs(chopper.speed);
  const unsigned opening_count = chopper.edge_count / 2;
  const double shortest_path = chopper.path;
  const double longest_path = chopper.path + path_spread;

  /* A spread wide enough to reach from one turn into the next stops the pieces being
   * disjoint, and then areas double-count. Test it against the slowest neutron present. */
  if (path_spread > 0.0 && set->count) {
    double inverse_velocity_maximum = 0.0;
    for (unsigned i = 0; i < set->count; ++i) {
      double high = 0.0;
      chopper_polygon_extent(&set->polygon[i], 1.0, 0.0, NULL, &high);
      if (high > inverse_velocity_maximum) inverse_velocity_maximum = high;
    }
    for (unsigned w = 0; w < opening_count; ++w) {
      double lo = 0, hi = 0;
      chopper_opening_times(chopper, w, &lo, &hi);
      if (hi - lo >= tau) continue;  /* never shut, so nothing to overlap into */
      if (path_spread * inverse_velocity_maximum >= tau - (hi - lo)) {
        printf("A path spread of %g m reaches from one turn of this chopper into the next "
               "for an inverse velocity of %g s/m; the transmitted pieces would overlap "
               "and their areas would count twice. Narrow the inverse velocity range or "
               "the spread.\n", path_spread, inverse_velocity_maximum);
        chopper_polygon_set_free(set);
        return 0;
      }
    }
  }

  chopper_polygon_set out = chopper_polygon_set_empty();
  for (unsigned i = 0; i < set->count; ++i) {
    /* The earliest and latest this polygon can reach the disk: the near path bounds the
     * first and the far path the second. */
    double u_min = 0.0, u_max = 0.0, ignored = 0.0;
    chopper_polygon_extent(&set->polygon[i], shortest_path, 1.0, &u_min, &ignored);
    chopper_polygon_extent(&set->polygon[i], longest_path, 1.0, &ignored, &u_max);

    for (unsigned w = 0; w < opening_count; ++w) {
      double lo = 0, hi = 0;
      chopper_opening_times(chopper, w, &lo, &hi);
      /* Each failure below unwinds by hand rather than jumping to one exit. This file is
       * copied verbatim into a McStas instrument beside every other library the
       * instrument includes, and while a label cannot collide with another function's
       * label -- labels have function scope -- it can be eaten by a macro. Anything that
       * defines `failed` would rewrite both the label and the jumps to it, which is the
       * same way `PI` bites a caller that defines its own. Three lines, three times. */
      if (hi - lo >= tau) {
        /* Open for at least a whole turn: it constrains nothing, like one parked open. */
        if (!chopper_polygon_set_add(&out, &set->polygon[i])) {
          chopper_polygon_set_free(&out);
          chopper_polygon_set_free(set);
          return 0;
        }
        continue;
      }
      const long first = (long) floor((u_min - hi) / tau);
      const long last = (long) ceil((u_max - lo) / tau);
      for (long n = first; n <= last; ++n) {
        chopper_polygon piece = set->polygon[i];
        if (!chopper_polygon_clip_wedge(&piece, shortest_path, longest_path, 1.0,
                                        lo + (double) n * tau, hi + (double) n * tau)) {
          printf("A transmitted polygon needs more than %d vertices; raise "
                 "CHOPPER_POLYGON_MAX_VERTICES\n", CHOPPER_POLYGON_MAX_VERTICES);
          chopper_polygon_set_free(&out);
          chopper_polygon_set_free(set);
          return 0;
        }
        if (!chopper_polygon_set_add(&out, &piece)) {
          chopper_polygon_set_free(&out);
          chopper_polygon_set_free(set);
          return 0;
        }
      }
    }
  }
  chopper_polygon_set_free(set);
  *set = out;
  return 1;
}

int chopper_polygon_set_transmit_train(chopper_polygon_set * set, const unsigned count,
                                       const chopper_parameters * choppers,
                                       const double * path_spreads) {
  if (set == NULL || (count && choppers == NULL)) return 0;
  for (unsigned i = 0; i < count && set->count; ++i) {
    if (!chopper_polygon_set_transmit(set, choppers[i],
                                      path_spreads ? path_spreads[i] : 0.0)) return 0;
  }
  return 1;
}

range_set chopper_polygon_set_inverse_velocity_ranges(const chopper_polygon_set * set) {
  range_set ranges;
  ranges.count = 0;
  ranges.ranges = NULL;
  if (set == NULL || set->count == 0) return ranges;
  ranges.ranges = (range *) calloc(set->count, sizeof(range));
  if (ranges.ranges == NULL) {
    printf("Out of memory projecting %u polygons\n", set->count);
    return ranges;
  }
  for (unsigned i = 0; i < set->count; ++i) {
    chopper_polygon_extent(&set->polygon[i], 1.0, 0.0,
                           &ranges.ranges[i].minimum, &ranges.ranges[i].maximum);
  }
  ranges.count = set->count;
  /* Sorts in place and merges what overlaps or touches; may return the input. */
  const range_set merged = range_set_sort(ranges);
  if (merged.ranges != ranges.ranges && ranges.ranges) free(ranges.ranges);
  return merged;
}

#pragma acc routine seq
int chopper_polygon_contains(const chopper_polygon * polygon, const double inverse_velocity,
                             const double time) {
  if (polygon == NULL || polygon->count < 3) return 0;
  /* Convex, so an inside point is on the same side of every edge. Collecting both signs
   * rather than comparing against the first lets a collinear edge -- cross exactly zero,
   * which a point on the boundary gives -- count as neither. */
  int positive = 0, negative = 0;
  for (unsigned i = 0; i < polygon->count; ++i) {
    const chopper_point p = polygon->vertex[i];
    const chopper_point q = polygon->vertex[(i + 1) % polygon->count];
    const double cross = (q.inverse_velocity - p.inverse_velocity) * (time - p.time)
                       - (q.time - p.time) * (inverse_velocity - p.inverse_velocity);
    if (cross > 0.0) positive = 1;
    if (cross < 0.0) negative = 1;
    if (positive && negative) return 0;
  }
  return 1;
}

#pragma acc routine seq
int chopper_polygon_set_contains(const chopper_polygon_set * set,
                                 const double inverse_velocity, const double time) {
  if (set == NULL) return 0;
  for (unsigned i = 0; i < set->count; ++i) {
    if (chopper_polygon_contains(&set->polygon[i], inverse_velocity, time)) return 1;
  }
  return 0;
}

void chopper_polygon_sampler_empty(chopper_polygon_sampler * sampler) {
  if (sampler == NULL) return;
  sampler->count = 0;
  sampler->acceptance = 0.0;
  sampler->cumulative = NULL;
  sampler->origin = NULL;
  sampler->edge_a = NULL;
  sampler->edge_b = NULL;
}

void chopper_polygon_sampler_free(chopper_polygon_sampler * sampler) {
  if (sampler == NULL) return;
  if (sampler->cumulative) free(sampler->cumulative);
  if (sampler->origin) free(sampler->origin);
  if (sampler->edge_a) free(sampler->edge_a);
  if (sampler->edge_b) free(sampler->edge_b);
  chopper_polygon_sampler_empty(sampler);
}

chopper_polygon_sampler chopper_polygon_sampler_make(const chopper_polygon_set * set,
                                                     const double sampled_area) {
  chopper_polygon_sampler sampler;
  chopper_polygon_sampler_empty(&sampler);
  if (set == NULL || set->count == 0 || !(sampled_area > 0.0)) return sampler;

  unsigned triangles = 0;
  for (unsigned i = 0; i < set->count; ++i) {
    if (set->polygon[i].count >= 3) triangles += set->polygon[i].count - 2;
  }
  if (triangles == 0) return sampler;

  sampler.cumulative = (double *) calloc(triangles, sizeof(double));
  sampler.origin = (chopper_point *) calloc(triangles, sizeof(chopper_point));
  sampler.edge_a = (chopper_point *) calloc(triangles, sizeof(chopper_point));
  sampler.edge_b = (chopper_point *) calloc(triangles, sizeof(chopper_point));
  if (sampler.cumulative == NULL || sampler.origin == NULL
      || sampler.edge_a == NULL || sampler.edge_b == NULL) {
    printf("Out of memory building a sampler over %u triangles\n", triangles);
    chopper_polygon_sampler_free(&sampler);
    return sampler;
  }

  /* Fan each polygon from its first vertex. Convex, so every fan triangle is inside it. */
  double total = 0.0;
  unsigned t = 0;
  for (unsigned i = 0; i < set->count; ++i) {
    const chopper_polygon * polygon = &set->polygon[i];
    if (polygon->count < 3) continue;
    const chopper_point origin = polygon->vertex[0];
    for (unsigned v = 1; v + 1 < polygon->count; ++v) {
      const chopper_point a = polygon->vertex[v], b = polygon->vertex[v + 1];
      sampler.origin[t] = origin;
      sampler.edge_a[t].inverse_velocity = a.inverse_velocity - origin.inverse_velocity;
      sampler.edge_a[t].time = a.time - origin.time;
      sampler.edge_b[t].inverse_velocity = b.inverse_velocity - origin.inverse_velocity;
      sampler.edge_b[t].time = b.time - origin.time;
      total += fabs(sampler.edge_a[t].inverse_velocity * sampler.edge_b[t].time
                    - sampler.edge_b[t].inverse_velocity * sampler.edge_a[t].time) / 2.0;
      sampler.cumulative[t] = total;
      ++t;
    }
  }
  if (!(total > 0.0)) {
    chopper_polygon_sampler_free(&sampler);
    return sampler;
  }
  sampler.count = t;
  sampler.acceptance = total / sampled_area;
  for (unsigned i = 0; i < sampler.count; ++i) sampler.cumulative[i] /= total;
  /* Pin the last, rather than leave a draw above it to fall off the end of the search. */
  sampler.cumulative[sampler.count - 1] = 1.0;
  return sampler;
}

#pragma acc routine seq
void chopper_polygon_sampler_draw(const chopper_polygon_sampler * sampler,
                                  const double triangle_deviate, const double first_deviate,
                                  const double second_deviate,
                                  double * inverse_velocity, double * time) {
  if (sampler == NULL || sampler->count == 0) return;
  unsigned low = 0, high = sampler->count - 1;
  while (low < high) {
    const unsigned mid = low + (high - low) / 2;
    if (sampler->cumulative[mid] <= triangle_deviate) low = mid + 1; else high = mid;
  }
  /* Two uniform deviates land in the unit square; folding the far half back across the
   * diagonal puts them uniformly in the unit triangle, and the affine map carries that
   * to the real one. */
  double u = first_deviate, v = second_deviate;
  if (u + v > 1.0) { u = 1.0 - u; v = 1.0 - v; }
  if (inverse_velocity != NULL) {
    *inverse_velocity = sampler->origin[low].inverse_velocity
                      + u * sampler->edge_a[low].inverse_velocity
                      + v * sampler->edge_b[low].inverse_velocity;
  }
  if (time != NULL) {
    *time = sampler->origin[low].time + u * sampler->edge_a[low].time
          + v * sampler->edge_b[low].time;
  }
}

static void chopper_write_axes_to_file(FILE * file,
    const double * inverse_velocities, const unsigned inverse_velocity_count,
    const double * times, const unsigned time_count) {
  if (inverse_velocities && inverse_velocity_count) {
    fprintf(file, "# Inverse velocity edges (s/m):\n#");
    for (unsigned vi = 0; vi < inverse_velocity_count; ++vi) {
      fprintf(file, "%g ", inverse_velocities[vi]);
    }
    fprintf(file, "\n");
  }
  if (times && time_count) {
    fprintf(file, "# Time edges (s):\n#");
    for (unsigned ti = 0; ti < time_count; ++ti) {
      fprintf(file, "%g ", times[ti]);
    }
    fprintf(file, "\n");
  }
}

static FILE * chopper_open_file_for_writing(
  const char * directory, const char * filename, const char * extension, const char * path_sep
  ){
  const unsigned dlen = directory ? strlen(directory) : 0;
  const unsigned plen = path_sep ? strlen(path_sep) : 0;
  const unsigned flen = filename ? strlen(filename) : 0;
  const unsigned elen = extension ? strlen(extension) : 0;
  char * filepath = calloc(dlen + plen + flen + elen + 1, sizeof(char));
  int dp = 0, fp = 0;
  if (directory && dlen) {
    sprintf(filepath, "%s", directory);
    dp = strcmp(directory + dlen - plen, path_sep);
  }
  if (filename && flen) {
    fp = strcmp(filename, path_sep);
  }
  if (dp && fp) {
    // need to add a path separator between directory and filename
    sprintf(filepath + strlen(filepath), "%s",path_sep);
  }
  if (filename) {
    sprintf(filepath + strlen(filepath), "%s", filename);
  }
  if (elen && filename && extension && strcmp(filename + flen - elen, extension) != 0) {
    // need to add the extension
    sprintf(filepath + strlen(filepath), "%s", extension);
  }
  /* Truncating, not appending. A caller may write the same grid more than once in a run --
   * McStas saves on SIGUSR2 and carries on, then saves again at the end -- and appending
   * leaves the second copy nose to tail with the first in one file, which reads as a single
   * grid of twice the rows. */
  FILE * file = fopen(filepath, "w");
  if (file == NULL) {
    printf("Could not open file %s for writing\n", filepath);
  }
  free(filepath);
  return file;
}

int chopper_write_mask_to_file(
  const char * directory, const char * filename, const char * extension, const char * path_sep,
  const int * mask, const unsigned inverse_velocity_count, const unsigned time_count,
  const double * inverse_velocities, const double * times
){
  FILE * file = chopper_open_file_for_writing(directory, filename, extension, path_sep);
  if (file == NULL) {
    return -1;
  }
  fprintf(file, "# Chopper mask file generated by chopper-lib\n");
  chopper_write_axes_to_file(file, inverse_velocities, inverse_velocity_count+1, times, time_count+1);
  fprintf(file, "# Mask (rows: time bins, columns: inverse velocity bins):\n");

  for (unsigned ti = 0; ti < time_count; ++ti) {
    for (unsigned vi = 0; vi < inverse_velocity_count; ++vi) {
      const unsigned index = ti * inverse_velocity_count + vi;
      fprintf(file, "%d ", mask[index]);
    }
    fprintf(file, "\n");
  }
  fclose(file);
  return 0;
}

int chopper_write_polygons_to_file(
  const char * directory, const char * filename, const char * extension, const char * path_sep,
  const chopper_polygon_set * set, const chopper_polygon * sampled
) {
  FILE * file = chopper_open_file_for_writing(directory, filename, extension, path_sep);
  if (file == NULL) return 0;

  const double transmitted = chopper_polygon_set_area(set);
  const double sampled_area = sampled ? chopper_polygon_area(sampled) : 0.0;

  /* %.17g round-trips a double exactly, which is the point of writing the region rather
   * than a picture of it: a reader gets the vertices the calculation actually used. */
  fprintf(file, "{\n");
  fprintf(file, "  \"chopper_lib_version\": \"%d.%d.%d\",\n", CHOPPER_LIB_VERSION_MAJOR,
          CHOPPER_LIB_VERSION_MINOR, CHOPPER_LIB_VERSION_PATCH);
  fprintf(file, "  \"inverse_velocity_unit\": \"s/m\",\n");
  fprintf(file, "  \"time_unit\": \"s\",\n");
  if (sampled && sampled->count >= 3 && sampled_area > 0.0) {
    double iv_low = 0, iv_high = 0, t_low = 0, t_high = 0;
    chopper_polygon_extent(sampled, 1.0, 0.0, &iv_low, &iv_high);
    chopper_polygon_extent(sampled, 0.0, 1.0, &t_low, &t_high);
    fprintf(file, "  \"sampled\": {\"inverse_velocity\": [%.17g, %.17g], "
                  "\"time\": [%.17g, %.17g], \"area\": %.17g},\n",
            iv_low, iv_high, t_low, t_high, sampled_area);
    fprintf(file, "  \"acceptance\": %.17g,\n", transmitted / sampled_area);
  } else {
    fprintf(file, "  \"sampled\": null,\n");
    fprintf(file, "  \"acceptance\": null,\n");
  }
  fprintf(file, "  \"transmitted_area\": %.17g,\n", transmitted);

  range_set bands = chopper_polygon_set_inverse_velocity_ranges(set);
  fprintf(file, "  \"inverse_velocity_bands\": [");
  for (unsigned i = 0; i < bands.count; ++i) {
    fprintf(file, "%s[%.17g, %.17g]", i ? ", " : "",
            bands.ranges[i].minimum, bands.ranges[i].maximum);
  }
  fprintf(file, "],\n");
  if (bands.ranges) free(bands.ranges);

  fprintf(file, "  \"polygons\": [\n");
  for (unsigned i = 0; set != NULL && i < set->count; ++i) {
    const chopper_polygon * polygon = &set->polygon[i];
    fprintf(file, "    {\"area\": %.17g, \"vertices\": [", chopper_polygon_area(polygon));
    for (unsigned v = 0; v < polygon->count; ++v) {
      fprintf(file, "%s[%.17g, %.17g]", v ? ", " : "",
              polygon->vertex[v].inverse_velocity, polygon->vertex[v].time);
    }
    fprintf(file, "]}%s\n", (i + 1 < set->count) ? "," : "");
  }
  fprintf(file, "  ]\n}\n");

  fclose(file);
  return 1;
}

int chopper_write_total_to_file(
  const char * directory, const char * filename, const char * extension, const char * path_sep,
  const double * total, const unsigned inverse_velocity_count, const unsigned time_count,
  const double * inverse_velocities, const double * times
) {
  FILE * file = chopper_open_file_for_writing(directory, filename, extension, path_sep);
  if (file == NULL) {
    return -1;
  }
  fprintf(file, "# Chopper mask file generated by chopper-lib\n");
  chopper_write_axes_to_file(file, inverse_velocities, inverse_velocity_count + 1, times, time_count + 1);
  fprintf(file, "# Total signal distribution (rows: time bins, columns: inverse velocity bins):\n");

  for (unsigned ti = 0; ti < time_count; ++ti) {
    for (unsigned vi = 0; vi < inverse_velocity_count; ++vi) {
      const unsigned index = ti * inverse_velocity_count + vi;
      fprintf(file, "%g ", total[index]);
    }
    fprintf(file, "\n");
  }
  fclose(file);
  return 0;
}

#ifdef __cplusplus
} // end EXTERN "C"
#endif
