/* Disk choppers with more than one opening.
 *
 * `multi_chopper_inverse_velocity_windows` and its wavelength and envelope wrappers say
 * the same thing as the single-opening functions, for a disk whose openings are listed
 * individually by angle rather than summarised as one width. Three things have to hold
 * for that to be true, and the tests below are grouped by which one they pin:
 *
 *   - a one-window disk has to reproduce the single-opening answer exactly;
 *   - every opening on the disk has to contribute, at every rotation;
 *   - an opening has to be placed by its angle, which means by the *signed* speed, since
 *     a disk turning the other way reaches the same angle at a different time.
 *
 * The last of these is invisible to `single_to_multi_chopper`, whose one window is
 * symmetric about zero and so unchanged by reversal. It only shows up on a disk whose
 * openings are asymmetric, which is exactly what these functions exist to describe.
 */
#include <stdlib.h>
#include <string.h>
#include "chopper-lib.h"
#include "test_util.h"

#if !defined(CHOPPER_LIB_VERSION) || CHOPPER_LIB_VERSION < 30000
#error "These tests place window angles with the signed speed, which needs chopper-lib 3.0.0 or newer"
#endif

static const double PATH = 10.0;     /* m from the source to the chopper */
static const double SPEED = 14.0;    /* Hz */
static const double DELAY = 0.02;    /* s, when the zero-angle point is on the beam */
static const double IV_MIN = 0.0005; /* s/m */
static const double IV_MAX = 0.01;   /* s/m */

/* chopper-lib.c defines these privately for its wavelength wrappers; a test of a unit
 * conversion has to name the units it expects, so repeat them here rather than infer. */
#define TEST_V2K 1.58825361e-3
#define TEST_K2V 629.622368
#define TEST_PI 3.14159265358979323846

/* How wide, in inverse velocity, an opening of `degrees` is on a disk at `speed`. */
static double window_half_width(double degrees, double speed, double path) {
  return degrees / 2.0 / 360.0 / fabs(speed) / path;
}

static double centre_of(range r) { return (r.minimum + r.maximum) / 2.0; }
static double half_of(range r) { return (r.maximum - r.minimum) / 2.0; }

/* ---------------------------------------------------------------------------------
 * A one-window disk is the single-opening chopper, described differently.
 * ------------------------------------------------------------------------------- */

/* `single_to_multi_chopper` is the bridge between the two families of functions, so the
 * two have to agree on everything it produces -- including a disk turning backwards,
 * where they place the window by opposite conventions and land in the same place only
 * because a single opening is symmetric about its own centre.
 */
static void test_a_one_window_disk_matches_the_single_opening_answer(void) {
  TEST("a one-window disk gives the single-opening answer exactly");
  const chopper_parameters singles[4] = {
    { SPEED, DELAY, 3.6, PATH},
    {-SPEED, DELAY, 3.6, PATH},        /* reversed: symmetric window, same answer */
    { SPEED, 3.5 / SPEED, 12.0, PATH}, /* a delay of several periods */
    { 71.0, -0.004, 0.7, 162.0},       /* a fast, narrow, distant chopper */
  };

  for (int i = 0; i < 4; ++i) {
    range_set want = chopper_inverse_velocity_windows(1, &singles[i], IV_MIN, IV_MAX, 0.0);
    multi_chopper_parameters multi = single_to_multi_chopper(singles[i]);
    range_set got = multi_chopper_inverse_velocity_windows(1, &multi, IV_MIN, IV_MAX, 0.0);

    CHECK(want.count > 0);   /* an empty answer would match trivially */
    CHECK_EQUAL_INT(got.count, want.count);
    if (got.count == want.count) {
      for (unsigned w = 0; w < want.count; ++w) {
        CHECK_CLOSE(got.ranges[w].minimum, want.ranges[w].minimum, 1e-15);
        CHECK_CLOSE(got.ranges[w].maximum, want.ranges[w].maximum, 1e-15);
      }
    }
    free(multi.windows);
    if (want.ranges) free(want.ranges);
    if (got.ranges) free(got.ranges);
  }
}

/* The same equivalence for the envelope, which is what most callers actually use. */
static void test_a_one_window_disk_matches_the_single_opening_limits(void) {
  TEST("the envelope of a one-window disk matches the single-opening envelope");
  const chopper_parameters single = {SPEED, DELAY, 3.6, PATH};
  multi_chopper_parameters multi = single_to_multi_chopper(single);

  double want_lo = 0.0, want_hi = 0.0, got_lo = 0.0, got_hi = 0.0;
  const unsigned want = chopper_inverse_velocity_limits(&want_lo, &want_hi, 1, &single, IV_MIN, IV_MAX, 0.0);
  const unsigned got = multi_chopper_inverse_velocity_limits(&got_lo, &got_hi, 1, &multi, IV_MIN, IV_MAX, 0.0);

  CHECK(want > 0);
  CHECK_EQUAL_INT(got, want);
  CHECK_CLOSE(got_lo, want_lo, 1e-15);
  CHECK_CLOSE(got_hi, want_hi, 1e-15);
  free(multi.windows);
}

/* ---------------------------------------------------------------------------------
 * Every opening contributes, at every rotation.
 * ------------------------------------------------------------------------------- */

/* Six openings 60 degrees apart admit six times as many neutrons as one.
 *
 * The openings are evenly spaced, so they recur at a sixth of the disk period and the
 * answer is a comb six times as fine as the one-window answer -- the substitution the
 * `chopper_parameters` documentation describes, arrived at the long way round.
 *
 * This is also the case that used to run off the end of its own allocation: the range
 * array was sized for one range per rotation, but the loop fills one per rotation *per
 * window*, so any disk with more than one opening wrote past the end of the buffer.
 */
static void test_every_opening_on_an_evenly_spaced_disk_admits_a_window(void) {
  TEST("each opening of an evenly spaced disk admits its own window");
  const unsigned n_windows = 6;
  const double width = 4.0;               /* degrees, per opening */
  chopper_window windows[6];
  for (unsigned i = 0; i < n_windows; ++i) {
    windows[i].min = 360.0 * (double) i / (double) n_windows - width / 2.0;
    windows[i].max = 360.0 * (double) i / (double) n_windows + width / 2.0;
  }
  const multi_chopper_parameters disk = {
    .speed = SPEED, .delay = DELAY, .window_count = n_windows, .windows = windows, .path = PATH
  };

  range_set got = multi_chopper_inverse_velocity_windows(1, &disk, IV_MIN, IV_MAX, 0.0);

  /* Openings are on the beam at DELAY + k / (6 * SPEED) for every integer k, so in
   * inverse velocity they sit at (DELAY + k / 84) / PATH. Of those, k = -1 through 6
   * fall inside [IV_MIN, IV_MAX] whole; k = -2 and k = 7 fall outside it entirely. */
  const double spacing = 1.0 / ((double) n_windows * SPEED) / PATH;
  const double half = window_half_width(width, SPEED, PATH);
  CHECK_EQUAL_INT(got.count, 8);
  if (got.count == 8) {
    for (unsigned w = 0; w < 8; ++w) {
      const double k = (double) w - 1.0;
      CHECK_CLOSE(centre_of(got.ranges[w]), DELAY / PATH + k * spacing, 1e-15);
      CHECK_CLOSE(half_of(got.ranges[w]), half, 1e-15);
    }
  }
  if (got.ranges) free(got.ranges);
}

/* Openings sit where their angles put them, not merely n to a turn.
 *
 * An evenly spaced disk cannot tell "six openings" from "one opening at six times the
 * speed", so it cannot catch a fix that only multiplies the window count. Two openings
 * a quarter turn apart can: the gaps they leave alternate between a quarter period and
 * three quarters of one.
 */
static void test_openings_are_placed_by_angle(void) {
  TEST("openings land where their angles put them, not evenly spaced");
  chopper_window windows[2] = {{.min = -1.0, .max = 1.0}, {.min = 89.0, .max = 91.0}};
  const multi_chopper_parameters disk = {
    .speed = SPEED, .delay = DELAY, .window_count = 2, .windows = windows, .path = PATH
  };

  range_set got = multi_chopper_inverse_velocity_windows(1, &disk, IV_MIN, IV_MAX, 0.0);

  /* Beam crossings at DELAY + n / SPEED and DELAY + (90/360) / SPEED + n / SPEED. Inside
   * an arrival time of PATH * IV_MAX = 0.1 s that is 0.02, 0.02 + 1/56 and 0.02 + 1/14. */
  const double expected[3] = {
    DELAY / PATH,
    (DELAY + 0.25 / SPEED) / PATH,
    (DELAY + 1.0 / SPEED) / PATH,
  };
  CHECK_EQUAL_INT(got.count, 3);
  if (got.count == 3) {
    for (unsigned w = 0; w < 3; ++w) CHECK_CLOSE(centre_of(got.ranges[w]), expected[w], 1e-15);
    /* the gaps really are uneven: a quarter turn, then three quarters */
    const double first_gap = centre_of(got.ranges[1]) - centre_of(got.ranges[0]);
    const double second_gap = centre_of(got.ranges[2]) - centre_of(got.ranges[1]);
    CHECK_CLOSE(second_gap / first_gap, 3.0, 1e-9);
  }
  if (got.ranges) free(got.ranges);
}

/* The order the openings are listed in, and which way round each one's edges are given,
 * are both presentation. The answer is a set of ranges either way.
 */
static void test_the_answer_does_not_depend_on_how_the_windows_are_listed(void) {
  TEST("listing the openings differently does not change the answer");
  chopper_window ordered[3] = {{-2.0, 2.0}, {58.0, 62.0}, {178.0, 182.0}};
  /* shuffled, and with two of the three given maximum-first */
  chopper_window jumbled[3] = {{182.0, 178.0}, {-2.0, 2.0}, {62.0, 58.0}};
  const multi_chopper_parameters a = {
    .speed = SPEED, .delay = DELAY, .window_count = 3, .windows = ordered, .path = PATH};
  const multi_chopper_parameters b = {
    .speed = SPEED, .delay = DELAY, .window_count = 3, .windows = jumbled, .path = PATH};

  range_set want = multi_chopper_inverse_velocity_windows(1, &a, IV_MIN, IV_MAX, 0.0);
  range_set got = multi_chopper_inverse_velocity_windows(1, &b, IV_MIN, IV_MAX, 0.0);

  CHECK(want.count > 0);
  CHECK_EQUAL_INT(got.count, want.count);
  if (got.count == want.count) {
    for (unsigned w = 0; w < want.count; ++w) {
      CHECK_CLOSE(got.ranges[w].minimum, want.ranges[w].minimum, 1e-15);
      CHECK_CLOSE(got.ranges[w].maximum, want.ranges[w].maximum, 1e-15);
    }
  }
  if (want.ranges) free(want.ranges);
  if (got.ranges) free(got.ranges);
}

/* ---------------------------------------------------------------------------------
 * Which way the disk turns.
 * ------------------------------------------------------------------------------- */

/* Reversing the disk reflects its openings about the delay.
 *
 * `delay` fixes when the zero-angle point is on the beam, whichever way the disk turns.
 * An opening ahead of that point in the direction of travel is reached after it; reverse
 * the disk and the same opening is reached before it. A single opening centred on zero
 * cannot show this -- reflecting it maps it onto itself -- so it takes an asymmetric one.
 */
static void test_reversing_the_disk_reflects_asymmetric_openings(void) {
  TEST("reversing the disk reflects its openings about the delay");
  chopper_window window = {.min = 0.0, .max = 10.0};  /* entirely on one side of zero */
  const multi_chopper_parameters forward = {
    .speed = SPEED, .delay = DELAY, .window_count = 1, .windows = &window, .path = PATH};
  const multi_chopper_parameters reverse = {
    .speed = -SPEED, .delay = DELAY, .window_count = 1, .windows = &window, .path = PATH};

  range_set f = multi_chopper_inverse_velocity_windows(1, &forward, IV_MIN, IV_MAX, 0.0);
  range_set r = multi_chopper_inverse_velocity_windows(1, &reverse, IV_MIN, IV_MAX, 0.0);

  /* forward: the zero-angle edge arrives at DELAY and the opening runs 10 degrees later */
  const double reference = DELAY / PATH;
  const double span = 10.0 / 360.0 / SPEED / PATH;
  CHECK(f.count > 0);
  CHECK_EQUAL_INT(r.count, f.count);
  if (f.count > 0 && r.count == f.count) {
    CHECK_CLOSE(f.ranges[0].minimum, reference, 1e-15);
    CHECK_CLOSE(f.ranges[0].maximum, reference + span, 1e-15);
    /* reversed, the same opening precedes the zero-angle point instead of following it */
    CHECK_CLOSE(r.ranges[0].minimum, reference - span, 1e-15);
    CHECK_CLOSE(r.ranges[0].maximum, reference, 1e-15);
  }
  if (f.ranges) free(f.ranges);
  if (r.ranges) free(r.ranges);
}

/* Where the mask puts its allowed bins, as contiguous runs of inverse velocity. */
static unsigned mask_runs(const multi_chopper_parameters * disk, range * runs, unsigned max_runs) {
  const unsigned bins = 200000;
  double * edges = calloc(bins + 1, sizeof(double));
  int * mask = calloc(bins, sizeof(int));
  for (unsigned i = 0; i <= bins; ++i) edges[i] = IV_MAX * (double) i / (double) bins;
  const double times[2] = {0.0, 1.0e-12};   /* one vanishingly short emission bin */

  multi_chopper_inverse_velocity_time_mask(
    mask, bins, 1, edges, bins + 1, times, 2, disk, 1, 0 /* no growing */);

  unsigned found = 0;
  double start = 0.0;
  int inside = 0;
  for (unsigned i = 0; i < bins; ++i) {
    if (mask[i] == CHOPPER_MASK_INCLUDED && !inside) { inside = 1; start = edges[i]; }
    else if (mask[i] != CHOPPER_MASK_INCLUDED && inside) {
      inside = 0;
      if (found < max_runs) { runs[found].minimum = start; runs[found].maximum = edges[i]; }
      ++found;
    }
  }
  if (inside && found < max_runs) { runs[found].minimum = start; runs[found].maximum = IV_MAX; ++found; }
  free(edges);
  free(mask);
  return found;
}

/* The mask and the window list describe the same chopper.
 *
 * They are separate implementations -- one intersects ranges, the other rejects bins --
 * so they can drift apart. They did: the mask placed window angles with the unsigned
 * speed, which agrees for an opening symmetric about zero and disagrees for any other
 * one, in the direction of mirroring the whole disk.
 */
static void test_the_mask_agrees_with_the_window_list(void) {
  TEST("the time mask admits what the window list says it should");
  chopper_window windows[2] = {{.min = 0.0, .max = 10.0}, {.min = 100.0, .max = 104.0}};
  const double speeds[2] = {SPEED, -SPEED};

  for (int s = 0; s < 2; ++s) {
    const multi_chopper_parameters disk = {
      .speed = speeds[s], .delay = DELAY, .window_count = 2, .windows = windows, .path = PATH};

    range_set want = multi_chopper_inverse_velocity_windows(1, &disk, 0.0, IV_MAX, 0.0);
    range runs[16];
    const unsigned found = mask_runs(&disk, runs, 16);

    CHECK(want.count > 0);
    CHECK_EQUAL_INT(found, want.count);
    if (found == want.count) {
      /* a bin survives if any of its edges reaches an opening, so the mask is generous
       * by up to a bin at each end; the grid is fine enough that this is the tolerance */
      const double tolerance = 3.0 * IV_MAX / 200000.0;
      for (unsigned w = 0; w < found && w < 16; ++w) {
        CHECK_CLOSE(centre_of(runs[w]), centre_of(want.ranges[w]), tolerance);
        CHECK_CLOSE(half_of(runs[w]), half_of(want.ranges[w]), tolerance);
      }
    }
    if (want.ranges) free(want.ranges);
  }
}

/* ---------------------------------------------------------------------------------
 * Trains, envelopes and wavelengths.
 * ------------------------------------------------------------------------------- */

/* Two openings offer twice as many chances to pass; a second disk takes half of them back.
 *
 * The two-window disk admits a window every half period. The single-window disk behind it
 * admits one every whole period, at inverse velocities the first disk also admits, so the
 * train passes every other one of the first disk's windows and none of its own extra.
 */
static void test_a_train_keeps_only_the_openings_every_disk_admits(void) {
  TEST("a train of multi-window disks passes only their common openings");
  chopper_window two[2] = {{.min = -2.0, .max = 2.0}, {.min = 178.0, .max = 182.0}};
  chopper_window one[1] = {{.min = -1.0, .max = 1.0}};   /* half as wide */
  const multi_chopper_parameters train[2] = {
    {.speed = SPEED, .delay = DELAY, .window_count = 2, .windows = two, .path = PATH},
    {.speed = SPEED, .delay = DELAY, .window_count = 1, .windows = one, .path = PATH},
  };

  range_set alone = multi_chopper_inverse_velocity_windows(1, &train[0], IV_MIN, IV_MAX, 0.0);
  range_set both = multi_chopper_inverse_velocity_windows(2, train, IV_MIN, IV_MAX, 0.0);

  /* the first disk on its own admits DELAY/PATH and a window every 1/(2*SPEED)/PATH
   * either side of it: three of them fall inside [IV_MIN, IV_MAX] */
  CHECK_EQUAL_INT(alone.count, 3);
  /* behind the second disk only every other one survives */
  CHECK_EQUAL_INT(both.count, 2);
  if (both.count == 2) {
    const double half = window_half_width(2.0, SPEED, PATH);  /* the narrower opening */
    CHECK_CLOSE(centre_of(both.ranges[0]), DELAY / PATH, 1e-15);
    CHECK_CLOSE(centre_of(both.ranges[1]), (DELAY + 1.0 / SPEED) / PATH, 1e-15);
    for (unsigned w = 0; w < 2; ++w) CHECK_CLOSE(half_of(both.ranges[w]), half, 1e-15);
  }
  if (alone.ranges) free(alone.ranges);
  if (both.ranges) free(both.ranges);
}

/* A stationary disk is not a closed one: the window functions step over it. */
static void test_a_stationary_disk_is_ignored(void) {
  TEST("a disk that is not turning does not block anything");
  chopper_window window = {.min = -2.0, .max = 2.0};
  const multi_chopper_parameters train[2] = {
    {.speed = SPEED, .delay = DELAY, .window_count = 1, .windows = &window, .path = PATH},
    {.speed = 0.0, .delay = DELAY, .window_count = 1, .windows = &window, .path = PATH},
  };

  range_set turning = multi_chopper_inverse_velocity_windows(1, train, IV_MIN, IV_MAX, 0.0);
  range_set with_still = multi_chopper_inverse_velocity_windows(2, train, IV_MIN, IV_MAX, 0.0);

  CHECK(turning.count > 0);
  CHECK_EQUAL_INT(with_still.count, turning.count);
  if (with_still.count == turning.count) {
    for (unsigned w = 0; w < turning.count; ++w) {
      CHECK_CLOSE(with_still.ranges[w].minimum, turning.ranges[w].minimum, 1e-15);
      CHECK_CLOSE(with_still.ranges[w].maximum, turning.ranges[w].maximum, 1e-15);
    }
  }
  if (turning.ranges) free(turning.ranges);
  if (with_still.ranges) free(with_still.ranges);
}

/* The envelope spans the outermost windows, gaps and all.
 *
 * A multi-window disk makes this the common case rather than the exception, which is why
 * the count comes back alongside the bounds: two windows and a gap between them report
 * the same pair of numbers as one window covering the whole span.
 */
static void test_the_envelope_spans_the_outermost_windows(void) {
  TEST("the reported envelope spans every window, including the gaps");
  chopper_window windows[2] = {{.min = -2.0, .max = 2.0}, {.min = 178.0, .max = 182.0}};
  const multi_chopper_parameters disk = {
    .speed = SPEED, .delay = DELAY, .window_count = 2, .windows = windows, .path = PATH};

  range_set expected = multi_chopper_inverse_velocity_windows(1, &disk, IV_MIN, IV_MAX, 0.0);
  double lower = 0.0, upper = 0.0;
  const unsigned count = multi_chopper_inverse_velocity_limits(
    &lower, &upper, 1, &disk, IV_MIN, IV_MAX, 0.0);

  CHECK(expected.count > 1);   /* or there would be no gap to envelope */
  CHECK_EQUAL_INT(count, expected.count);
  if (expected.count) {
    CHECK_CLOSE(lower, expected.ranges[0].minimum, 1e-15);
    CHECK_CLOSE(upper, expected.ranges[expected.count - 1].maximum, 1e-15);
    /* the envelope is wider than the windows it covers, which is what the count warns of */
    CHECK(upper - lower > expected.ranges[0].maximum - expected.ranges[0].minimum);
  }
  if (expected.ranges) free(expected.ranges);
}

/* A chopper train that admits nothing reports nothing, and leaves the outputs alone. */
static void test_a_train_that_admits_nothing_reports_nothing(void) {
  TEST("a train that passes nothing returns no windows and writes no bounds");
  chopper_window a = {.min = -0.5, .max = 0.5};
  chopper_window b = {.min = -0.5, .max = 0.5};
  /* the second disk is a quarter period out of step with the first, at the same distance */
  const multi_chopper_parameters train[2] = {
    {.speed = SPEED, .delay = DELAY, .window_count = 1, .windows = &a, .path = PATH},
    {.speed = SPEED, .delay = DELAY + 0.25 / SPEED, .window_count = 1, .windows = &b, .path = PATH},
  };

  double lower = -1.0, upper = -1.0;
  const unsigned count = multi_chopper_inverse_velocity_limits(
    &lower, &upper, 2, train, IV_MIN, IV_MAX, 0.0);

  CHECK_EQUAL_INT(count, 0);
  /* the documented contract: the bounds are only set when the count is non-zero */
  CHECK_CLOSE(lower, -1.0, 0.0);
  CHECK_CLOSE(upper, -1.0, 0.0);
}

/* The wavelength wrapper is the inverse-velocity answer in different units. */
static void test_wavelength_limits_are_the_inverse_velocity_limits_converted(void) {
  TEST("the wavelength envelope is the inverse velocity envelope, converted");
  chopper_window window = {.min = -2.0, .max = 2.0};
  const multi_chopper_parameters disk = {
    .speed = SPEED, .delay = DELAY, .window_count = 1, .windows = &window, .path = PATH};
  const double lambda_min = 1.0, lambda_max = 20.0;   /* angstrom */

  double lo = 0.0, hi = 0.0;
  const unsigned count = multi_chopper_wavelength_limits(
    &lo, &hi, 1, &disk, lambda_min, lambda_max, 0.0);

  /* Openings are at (DELAY + n / SPEED) / PATH in inverse velocity; over 1 to 20 A --
   * 2.53e-4 to 5.06e-3 s/m -- only n = 0 is in range, so the envelope is that one window. */
  const double centre = DELAY / PATH;
  const double half = window_half_width(4.0, SPEED, PATH);
  CHECK_EQUAL_INT(count, 1);
  CHECK_CLOSE(lo, (centre - half) * TEST_K2V * 2 * TEST_PI, 1e-9);
  CHECK_CLOSE(hi, (centre + half) * TEST_K2V * 2 * TEST_PI, 1e-9);
  /* 0.002 s/m is 500 m/s, which is a shade under 8 A -- a sanity check on the constants */
  CHECK_CLOSE((lo + hi) / 2.0, 7.9125, 1e-3);

  /* and the same chopper asked in inverse velocity gives the bounds those came from */
  double iv_lo = 0.0, iv_hi = 0.0;
  const unsigned iv_count = multi_chopper_inverse_velocity_limits(
    &iv_lo, &iv_hi, 1, &disk,
    lambda_min * TEST_V2K / 2 / TEST_PI, lambda_max * TEST_V2K / 2 / TEST_PI, 0.0);
  CHECK_EQUAL_INT(iv_count, count);
  CHECK_CLOSE(lo, iv_lo * TEST_K2V * 2 * TEST_PI, 1e-12);
  CHECK_CLOSE(hi, iv_hi * TEST_K2V * 2 * TEST_PI, 1e-12);
}

/* A pulse of finite length stretches every window towards shorter inverse velocity.
 *
 * A neutron emitted at the end of the pulse has to fly faster to reach the same opening
 * at the same time, so the lower edge of each window moves down by the pulse length over
 * the flight path while the upper edge, which belongs to a neutron emitted at t = 0,
 * stays put.
 */
static void test_a_longer_pulse_widens_every_window_downwards(void) {
  TEST("a finite emission time widens each window towards shorter inverse velocity");
  chopper_window windows[2] = {{.min = -2.0, .max = 2.0}, {.min = 118.0, .max = 122.0}};
  const multi_chopper_parameters disk = {
    .speed = SPEED, .delay = DELAY, .window_count = 2, .windows = windows, .path = PATH};
  const double emission = 2.86e-3;   /* s, roughly an ESS pulse */

  range_set sharp = multi_chopper_inverse_velocity_windows(1, &disk, IV_MIN, IV_MAX, 0.0);
  range_set spread = multi_chopper_inverse_velocity_windows(1, &disk, IV_MIN, IV_MAX, emission);

  CHECK(sharp.count > 0);
  CHECK_EQUAL_INT(spread.count, sharp.count);
  if (spread.count == sharp.count) {
    for (unsigned w = 0; w < sharp.count; ++w) {
      CHECK_CLOSE(spread.ranges[w].maximum, sharp.ranges[w].maximum, 1e-15);
      CHECK_CLOSE(spread.ranges[w].minimum,
                  sharp.ranges[w].minimum - emission / PATH, 1e-15);
    }
  }
  if (sharp.ranges) free(sharp.ranges);
  if (spread.ranges) free(spread.ranges);
}

/* The single-opening wavelength wrapper and the multi-opening one are the same wrapper.
 *
 * Both divide the work the same way -- convert to inverse velocity, intersect, convert
 * back -- so a one-window disk has to come out of either in the same place. This also
 * pins the conversion itself: an error in V2K or K2V would move both, but an error in
 * only one of the two wrappers moves one.
 */
static void test_the_two_wavelength_wrappers_agree(void) {
  TEST("both wavelength wrappers give the same envelope for the same chopper");
  const chopper_parameters single = {SPEED, DELAY, 4.0, PATH};
  multi_chopper_parameters multi = single_to_multi_chopper(single);
  const double lambda_min = 1.0, lambda_max = 20.0;

  double want_lo = 0.0, want_hi = 0.0, got_lo = 0.0, got_hi = 0.0;
  const unsigned want = chopper_wavelength_limits(
    &want_lo, &want_hi, 1, &single, lambda_min, lambda_max, 0.0);
  const unsigned got = multi_chopper_wavelength_limits(
    &got_lo, &got_hi, 1, &multi, lambda_min, lambda_max, 0.0);

  CHECK(want > 0);
  CHECK_EQUAL_INT(got, want);
  CHECK_CLOSE(got_lo, want_lo, 1e-12);
  CHECK_CLOSE(got_hi, want_hi, 1e-12);
  /* and the answer is in angstrom, not seconds per metre */
  CHECK(want_lo > 1.0 && want_hi < 20.0);
  free(multi.windows);
}

/* A disk with no openings in it is a beam stop.
 *
 * It is a degenerate description rather than a real chopper, but the two implementations
 * have to make the same of it or a caller that builds windows dynamically gets one answer
 * from the envelope and the opposite one from the mask.
 */
static void test_a_disk_with_no_openings_passes_nothing(void) {
  TEST("a disk with no openings blocks everything, both ways of asking");
  const multi_chopper_parameters blocked = {
    .speed = SPEED, .delay = DELAY, .window_count = 0, .windows = NULL, .path = PATH};

  range_set windows = multi_chopper_inverse_velocity_windows(1, &blocked, IV_MIN, IV_MAX, 0.0);
  CHECK_EQUAL_INT(windows.count, 0);
  if (windows.ranges) free(windows.ranges);

  range runs[4];
  CHECK_EQUAL_INT(mask_runs(&blocked, runs, 4), 0);
}

int main(void) {
  test_a_one_window_disk_matches_the_single_opening_answer();
  test_a_one_window_disk_matches_the_single_opening_limits();
  test_every_opening_on_an_evenly_spaced_disk_admits_a_window();
  test_openings_are_placed_by_angle();
  test_the_answer_does_not_depend_on_how_the_windows_are_listed();
  test_reversing_the_disk_reflects_asymmetric_openings();
  test_the_mask_agrees_with_the_window_list();
  test_a_train_keeps_only_the_openings_every_disk_admits();
  test_a_stationary_disk_is_ignored();
  test_the_envelope_spans_the_outermost_windows();
  test_a_train_that_admits_nothing_reports_nothing();
  test_wavelength_limits_are_the_inverse_velocity_limits_converted();
  test_the_two_wavelength_wrappers_agree();
  test_a_disk_with_no_openings_passes_nothing();
  test_a_longer_pulse_widens_every_window_downwards();
  return chopper_test_report();
}
