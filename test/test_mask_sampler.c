/* What `chopper_mask_sampler_make` and `_draw` do with a mask whose answer is known.
 *
 * The sampler exists so a source can draw straight from the cells a chopper train allows
 * instead of drawing everywhere and throwing most of it away. That is only the same
 * measurement if two things hold, and both are checked here:
 *
 *   - `acceptance` is the allowed area over the sampled area. It is the factor the ray
 *     weight is multiplied by, so an error in it is a scale error on every monitor in the
 *     instrument and nothing else would catch it.
 *   - the draws are uniform over the allowed cells and land nowhere else. Not merely
 *     inside the mask -- *uniform*, because a source that favours one end of its allowed
 *     band reports the wrong spectrum with the right total.
 *
 * The grids here are deliberately not commensurate with the sampled region: the last row
 * and the last column stick out past it, which is what `ceil`-sized grids do in
 * `Masked_ESS_butterfly` and what the clipping in `_make` is for. A sampler that weighted
 * those cells whole would pass every other check in this file.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "chopper-lib.h"
#include "test_util.h"

/* A deterministic uniform generator, so a failure here is reproducible. Any decent LCG
 * does; this is the one from Numerical Recipes' `ranqd1`, taken as a fraction of 2^32. */
static unsigned long long lcg_state = 88172645463325252ull;
static double uniform(void) {
  lcg_state = 1664525ull * lcg_state + 1013904223ull;
  return (double)((lcg_state >> 16) & 0xffffffffull) / 4294967296.0;
}

/* Uniform edges, `count + 1` of them, starting at `low` and stepping by `step`. */
static double * edges(const double low, const double step, const unsigned count) {
  double * e = calloc(count + 1, sizeof(double));
  for (unsigned i = 0; i <= count; ++i) e[i] = low + (double)i * step;
  return e;
}

int main(void) {
  /* A 4 x 3 mask (4 inverse velocity bins, 3 time bins) with a known allowed set. Bins are
   * 1 s/m and 1 s wide from zero, but the sampled region is only 3.5 by 2.5 -- so the
   * fourth column and the third row are half outside it. */
  const unsigned nv = 4, nt = 3;
  double * ive = edges(0.0, 1.0, nv);
  double * te = edges(0.0, 1.0, nt);
  const double v_min = 0.0, v_range = 3.5;
  const double t_min = 0.0, t_range = 2.5;

  int mask[12];
  for (unsigned i = 0; i < nv * nt; ++i) mask[i] = CHOPPER_MASK_EXCLUDED;
  /* index is ti * nv + vi, as chopper_inverse_velocity_time_mask leaves it */
  mask[0 * nv + 1] = CHOPPER_MASK_INCLUDED;  /* whole:      1 x 1   = 1.0  */
  mask[0 * nv + 3] = CHOPPER_MASK_INCLUDED;  /* clipped v:  0.5 x 1 = 0.5  */
  mask[2 * nv + 1] = CHOPPER_MASK_INCLUDED;  /* clipped t:  1 x 0.5 = 0.5  */
  mask[2 * nv + 3] = CHOPPER_MASK_INCLUDED;  /* both:     0.5 x 0.5 = 0.25 */
  const double allowed_area = 1.0 + 0.5 + 0.5 + 0.25;

  chopper_mask_sampler sampler = chopper_mask_sampler_make(
    mask, nv, nt, ive, te, v_min, v_range, t_min, t_range);

  TEST("a clipped mask reports the area it really covers");
  CHECK_EQUAL_INT(sampler.count, 4);
  CHECK_CLOSE(sampler.acceptance, allowed_area / (v_range * t_range), 1e-12);
  /* The whole-cell count would give 4 / 12; the clipped area is 2.25 / 8.75. A sampler
   * that skipped the clipping lands on the first, so pin that they differ. */
  CHECK(fabs(sampler.acceptance - 4.0 / 12.0) > 1e-3);

  TEST("the cell distribution is normalised and increasing");
  for (unsigned i = 1; i < sampler.count; ++i) CHECK(sampler.cumulative[i] > sampler.cumulative[i - 1]);
  CHECK_CLOSE(sampler.cumulative[sampler.count - 1], 1.0, 0.0);

  TEST("every cell keeps only the part inside the sampled region");
  for (unsigned i = 0; i < sampler.count; ++i) {
    CHECK(sampler.inverse_velocity_low[i] >= v_min);
    CHECK(sampler.inverse_velocity_low[i] + sampler.inverse_velocity_width[i] <= v_min + v_range + 1e-12);
    CHECK(sampler.time_low[i] >= t_min);
    CHECK(sampler.time_low[i] + sampler.time_width[i] <= t_min + t_range + 1e-12);
    CHECK(sampler.inverse_velocity_width[i] > 0.0);
    CHECK(sampler.time_width[i] > 0.0);
  }

  /* Draw a great many points and bin them back onto the mask grid. Every one has to land
   * in an allowed cell, and the share landing in each has to be that cell's share of the
   * allowed area -- which is the whole claim the weight correction rests on. */
  const unsigned draws = 400000;
  double landed[12];
  for (unsigned i = 0; i < nv * nt; ++i) landed[i] = 0.0;
  unsigned outside = 0;
  for (unsigned d = 0; d < draws; ++d) {
    double iv = -1.0, t = -1.0;
    chopper_mask_sampler_draw(&sampler, uniform(), uniform(), uniform(), &iv, &t);
    const int vi = (int)floor((iv - v_min) / 1.0);
    const int ti = (int)floor((t - t_min) / 1.0);
    if (vi < 0 || vi >= (int)nv || ti < 0 || ti >= (int)nt
        || iv > v_min + v_range || t > t_min + t_range) {
      ++outside;
      continue;
    }
    landed[ti * nv + vi] += 1.0;
  }

  TEST("no draw lands outside the sampled region");
  CHECK_EQUAL_INT(outside, 0);

  TEST("no draw lands in a masked-out cell");
  for (unsigned i = 0; i < nv * nt; ++i) {
    if (mask[i] == CHOPPER_MASK_EXCLUDED) CHECK_EQUAL_INT((long)landed[i], 0);
  }

  TEST("draws fall in proportion to the clipped area of each cell");
  /* Poisson on ~ n p; five sigma on the smallest share is a couple of hundred counts out
   * of four hundred thousand, so this is tight without being flaky. */
  const double expected[12] = {
    0, 1.0, 0, 0.5,
    0, 0,   0, 0,
    0, 0.5, 0, 0.25
  };
  for (unsigned i = 0; i < nv * nt; ++i) {
    if (mask[i] == CHOPPER_MASK_EXCLUDED) continue;
    const double want = draws * expected[i] / allowed_area;
    CHECK_CLOSE(landed[i], want, 5.0 * sqrt(want));
  }

  TEST("within one cell the draw is flat, not merely inside it");
  /* The widest whole cell, split in half: an implementation that placed the point at a
   * cell edge, or scaled the deviate by the unclipped width, fails here and nowhere else. */
  unsigned lower = 0, upper = 0;
  for (unsigned d = 0; d < draws; ++d) {
    double iv = -1.0, t = -1.0;
    chopper_mask_sampler_draw(&sampler, 0.0, uniform(), uniform(), &iv, &t);
    if (iv < sampler.inverse_velocity_low[0] + sampler.inverse_velocity_width[0] / 2.0) ++lower;
    else ++upper;
  }
  CHECK_CLOSE((double)lower, draws / 2.0, 5.0 * sqrt(draws / 4.0));
  CHECK_EQUAL_INT(lower + upper, draws);

  chopper_mask_sampler_free(&sampler);

  TEST("freeing empties the sampler, and freeing twice is not an error");
  CHECK_EQUAL_INT(sampler.count, 0);
  CHECK(sampler.cumulative == NULL);
  chopper_mask_sampler_free(&sampler);
  chopper_mask_sampler_free(NULL);

  TEST("a mask that allows nothing yields a sampler that offers nothing");
  int shut[12];
  for (unsigned i = 0; i < nv * nt; ++i) shut[i] = CHOPPER_MASK_EXCLUDED;
  chopper_mask_sampler none = chopper_mask_sampler_make(
    shut, nv, nt, ive, te, v_min, v_range, t_min, t_range);
  CHECK_EQUAL_INT(none.count, 0);
  CHECK_CLOSE(none.acceptance, 0.0, 0.0);
  /* Drawing from it leaves the outputs alone rather than reading off the end of nothing */
  double untouched_iv = -7.0, untouched_t = -7.0;
  chopper_mask_sampler_draw(&none, 0.5, 0.5, 0.5, &untouched_iv, &untouched_t);
  CHECK_CLOSE(untouched_iv, -7.0, 0.0);
  CHECK_CLOSE(untouched_t, -7.0, 0.0);
  chopper_mask_sampler_free(&none);

  TEST("a mask allowing everything accepts everything");
  int open[12];
  for (unsigned i = 0; i < nv * nt; ++i) open[i] = CHOPPER_MASK_INCLUDED;
  chopper_mask_sampler all = chopper_mask_sampler_make(
    open, nv, nt, ive, te, v_min, v_range, t_min, t_range);
  CHECK_CLOSE(all.acceptance, 1.0, 1e-12);
  chopper_mask_sampler_free(&all);

  free(ive);
  free(te);
  return chopper_test_report();
}
