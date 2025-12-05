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
int compare_sorted_ranges(const void * ptr_a, const void * ptr_b){
  return compare_ranges((range *) ptr_a, (range *) ptr_b);
}
/******************************** range_set functions ************************************/
range_set range_set_sort(range_set s){
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

int range_intersects_ranges(const range r, const range_set rs){
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

/********************** local helper functions *************************/
// void print_help(const char * program_name) {
//   printf("Usage: %s [parameter name]=value ...\n", program_name);
//   printf("\tValid parameter names: xxNtype for xx in (ps, fo, bw), N in (1, 2), and type in (speed, phase)\n");
//   exit(0);
// }
//
// void print_one(const char * name, const double value) {
//   printf("  %s = %+10.4f\n", name, value);
// }
/*
 * toff=fabs(t-atan2(x,yprime)/omega - delay - (jitter ? jitter*randnorm():0));
   // does neutron hit outside slit? *
   if (fmod(toff+To/2.0,Tg)>To) ABSORB;
 */

/****************** chopper train functionality ************************/
range_set chopper_inverse_velocity_windows(const unsigned count, const chopper_parameters * choppers,
                                           const double inv_v_min, const double inv_v_max,
                                           const double latest_emission){
  range_set limits;
  limits.count = 1;
  limits.ranges = calloc(1, sizeof(range));
  limits.ranges[0].minimum = inv_v_min;
  limits.ranges[0].maximum = inv_v_max;
  for (unsigned i=0; i<count && limits.count; ++i) if (choppers[i].speed) {
    // the period of the chopper is a positive time
    const double tau = 1.0 / fabs(choppers[i].speed);
    // the delta time of the chopper is half the time it takes to rotate through the angle of the slit
    const double dt = choppers[i].angle / 360.0 / 2.0 * tau;
    // to match McStas DiskChopper, the delay time depends on the absolute value of the speed?
    const double t0 = choppers[i].phase / 360.0 / fabs(choppers[i].speed);
    // find the smallest n for which (t0 + dt + n * tau) / d >= inv_v_min
    const int n_min = (int) floor((choppers[i].path * inv_v_min - t0 - dt) / tau);
    // find the largest n for which (t0 - dt + n * tau) / d <= inv_v_max
    const int n_max = (int) ceil((choppers[i].path * inv_v_max - t0 + dt) / tau);
    // collect the ranges for each of the n in (n_min, n_max) into a set:
    range_set ith;
    ith.count = (unsigned)(n_max - n_min + 1);
    ith.ranges = calloc(ith.count, sizeof(range));
    if (ith.ranges == NULL) {
      printf("Out of memory\n");
      exit(-1);
    }
    unsigned c=0;
    for (unsigned j=0; j < ith.count; ++j) {
      const double n_tau = tau * (double) (n_min + (int) j);
      // let the minimum 1/v come from the *end* of the pulse:
      double j_min = (t0 - dt + n_tau - latest_emission) / choppers[i].path;
      double j_max = (t0 + dt + n_tau) / choppers[i].path;
      j_min = j_min < inv_v_min ? inv_v_min : j_min > inv_v_max ? inv_v_max : j_min;
      j_max = j_max < inv_v_min ? inv_v_min : j_max > inv_v_max ? inv_v_max : j_max;
      if (j_min < j_max) {
        ith.ranges[c].minimum = j_min;
        ith.ranges[c++].maximum = j_max;
      }
    }
    ith.count = c;
    // find the intersection of this chopper with the running set
    range_set new_limits = range_intersection(limits, ith);
    // clean-up allocated memory, ensuring limits or ith is not erased if transferred to new_limits:
    if ((!new_limits.count || new_limits.ranges != limits.ranges) && limits.ranges)  free(limits.ranges);
    if ((!new_limits.count || new_limits.ranges != ith.ranges) && ith.ranges) free(ith.ranges);
    // rename the new limits in preparation of returning or the next loop
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

// Use the McStas defines if possible, or define them ourselves
#ifndef V2K
#define V2K 1.58825361e-3     /* Convert v[m/s] to k[1/AA] */
#endif
#ifndef K2V
#define K2V 629.622368        /* Convert k[1/AA] to v[m/s] */
#endif
#ifndef PI
#define PI 3.14159265358979323846
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

multi_chopper_parameters single_to_multi_chopper(const chopper_parameters single) {
  chopper_window *  windows = calloc(1, sizeof(chopper_window)); // !leaked if not cleaned up externally
  windows->min = -single.angle/2.0;
  windows->max = single.angle/2.0;
  const multi_chopper_parameters multi = {single.speed, single.phase, 1, windows, single.path};
  return multi;
}

int_range multi_chopper_rotation_limits(const multi_chopper_parameters chopper, const range time_range) {
  int_range rotations = {1, -1};
  if (chopper.window_count < 1 || chopper.windows == NULL) {
    return rotations;
  }
  range windows_range = {chopper.windows[0].min, chopper.windows[0].max};
  for (unsigned i = 1; i<chopper.window_count; ++i) {
    if (chopper.windows[i].min < windows_range.minimum) windows_range.minimum = chopper.windows[i].min;
    if (chopper.windows[i].max > windows_range.maximum) windows_range.maximum = chopper.windows[i].max;
  }
  const double t0 = chopper.phase / 360.0 / fabs(chopper.speed);
  // Convert the windows full range from angle to fractional rotations
  windows_range.minimum = windows_range.minimum / 360.0; // rotations
  windows_range.maximum = windows_range.maximum / 360.0; // rotations
  // find the number of rotations needed to place the maximum angle _before_ the earliest time:
  rotations.minimum = (int) floor(chopper.speed * time_range.minimum - t0 - windows_range.maximum) - 1;
  // and the number of rotations needed to place the minimum angle _after_ the latest time:
  rotations.maximum = (int) ceil(chopper.speed * time_range.maximum - t0 - windows_range.minimum) + 1;
  return rotations;
}

unsigned multi_chopper_inverse_velocity_time_mask(
  int *mask, const unsigned mask_inverse_velocity_count, const unsigned mask_time_count,
  const double *inverse_velocities, const unsigned inverse_velocity_count,
  const double *times, const unsigned time_count,
  const multi_chopper_parameters *choppers, const unsigned chopper_count,
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
    const double tau = 1.0 / fabs(choppers[ci].speed);
    const double t0 = choppers[ci].phase * tau / 360.0;
    const range time_range = {times[0] + choppers[ci].path * inverse_velocities[0], times[time_count - 1] + choppers[ci].path * inverse_velocities[inverse_velocity_count - 1]};
    const int_range rotations = multi_chopper_rotation_limits(choppers[ci], time_range);
    if (rotations.maximum < rotations.minimum) continue; // No possible rotations
    // build the ranges of allowed times for this chopper
    range_set allowed_times;
    allowed_times.count = (rotations.maximum - rotations.minimum + 1) * choppers[ci].window_count;
    allowed_times.ranges = calloc(allowed_times.count, sizeof(range));
    unsigned c = 0;
    for (int n = rotations.minimum; n <= rotations.maximum; ++n) {
      for (unsigned w = 0; w < choppers[ci].window_count; ++w) {
        allowed_times.ranges[c].minimum = t0 + ((double) n + choppers[ci].windows[w].min / 360.0) * tau;
        allowed_times.ranges[c++].maximum = t0 + ((double) n + choppers[ci].windows[w].max / 360.0) * tau;
      }
    }
    // Now check each bin edge against the allowed times
    for (unsigned ti = 0; ti < time_count; ++ti) {
      for (unsigned vi = 0; vi < mask_inverse_velocity_count; ++ vi) {
        if (time_edges[ti * mask_inverse_velocity_count + vi] == 1) {
          const range edge= {
            times[ti] + choppers[ci].path * inverse_velocities[vi],
            times[ti] + choppers[ci].path * inverse_velocities[vi + 1]
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
            times[ti] + choppers[ci].path * inverse_velocities[vi],
            times[ti + 1] + choppers[ci].path * inverse_velocities[vi]
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


unsigned chopper_inverse_velocity_time_mask(
    int *mask, const unsigned mask_inverse_velocity_count, const unsigned mask_time_count,
    const double *inverse_velocities, const unsigned inverse_velocity_count,
    const double *times, const unsigned time_count,
    const chopper_parameters *choppers, const unsigned chopper_count,
    const int grow_mask
    ) {
  multi_chopper_parameters * multi_choppers = calloc(chopper_count, sizeof(multi_chopper_parameters));
  if (multi_choppers == NULL) {
    printf("Out of memory\n");
    exit(-1);
  }
  for (unsigned i = 0; i < chopper_count; ++i) {
    multi_choppers[i] = single_to_multi_chopper(choppers[i]);
  }
  const unsigned allowed_bins = multi_chopper_inverse_velocity_time_mask(
    mask, mask_inverse_velocity_count,  mask_time_count,
    inverse_velocities, inverse_velocity_count,
    times, time_count,
    multi_choppers, chopper_count,
    grow_mask
  );
  for (unsigned i = 0; i < chopper_count; ++i) {
    if (multi_choppers[i].windows) free(multi_choppers[i].windows);
  }
  free(multi_choppers);
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

void chopper_write_axes_to_file(FILE * file,
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

FILE * chopper_open_file_for_writing(
  const char * directory, const char * filename, const char * extension, const char * path_sep
  ){
  unsigned dlen = directory ? strlen(directory) : 0;
  unsigned plen = path_sep ? strlen(path_sep) : 0;
  unsigned flen = filename ? strlen(filename) : 0;
  unsigned elen = extension ? strlen(extension) : 0;
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
  FILE * file = fopen(filepath, "aw");
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
