/* What a chopper's `delay` means, and that the library agrees with it everywhere.
 *
 * A chopper opening is on the path at `delay`, and every period thereafter. That is a
 * statement about time alone: it does not mention how fast the disk turns, or which way.
 * Each test below pins one consequence of saying it that way, and each of them fails
 * against the `phase` formulation this library used before version 2.0.0.
 */
#include <stdlib.h>
#include "chopper-lib.h"
#include "test_util.h"

#if !defined(CHOPPER_LIB_VERSION) || CHOPPER_LIB_VERSION < 20000
#error "These tests describe delay-based choppers, which need chopper-lib 2.0.0 or newer"
#endif

static const double PATH = 10.0;    /* m from the source to the chopper */
static const double ANGLE = 3.6;    /* deg, a hundredth of a revolution */

/* Where the mask puts its allowed bins, as a mean inverse velocity, or -1 if none. */
static double allowed_centroid(const int * mask, unsigned iv_bins, unsigned t_bins,
                               const double * iv_edges, unsigned * allowed_out) {
  double weight = 0.0, total = 0.0;
  unsigned allowed = 0;
  for (unsigned ti = 0; ti < t_bins; ++ti) {
    for (unsigned vi = 0; vi < iv_bins; ++vi) {
      if (mask[ti * iv_bins + vi] == CHOPPER_MASK_INCLUDED) {
        const double centre = (iv_edges[vi] + iv_edges[vi + 1]) / 2.0;
        weight += centre;
        total += 1.0;
        ++allowed;
      }
    }
  }
  if (allowed_out) *allowed_out = allowed;
  return total > 0 ? weight / total : -1.0;
}

/* A single chopper against a fine inverse-velocity grid and one emission-time bin. */
static unsigned mask_one_chopper(chopper_parameters chopper, int * mask,
                                 double * iv_edges, unsigned iv_bins) {
  const double times[2] = {0.0, 1.0e-5};
  for (unsigned i = 0; i <= iv_bins; ++i) iv_edges[i] = 0.01 * (double) i / (double) iv_bins;
  return chopper_inverse_velocity_time_mask(
    mask, iv_bins, 1, iv_edges, iv_bins + 1, times, 2, &chopper, 1, 0 /* no growing */);
}

/* The window a delay names does not move when the disk spins faster.
 *
 * This is the whole reason for taking a delay rather than a phase: an angle has to be
 * divided by the speed to become a time, so the same phase at twice the speed named a
 * different window. The same delay must name the same window at any speed.
 */
static void test_the_delay_names_the_same_window_at_any_speed(void) {
  TEST("a delay names the same window however fast the disk turns");
  const double delay = 0.02;
  const double speeds[2] = {14.0, 28.0};

  for (int i = 0; i < 2; ++i) {
    const chopper_parameters chopper = {speeds[i], delay, ANGLE, PATH};
    range_set windows = chopper_inverse_velocity_windows(1, &chopper, 0.0005, 0.01, 0.0);

    /* the n=0 opening is centred on delay, so at delay/path in inverse velocity */
    const double expected_centre = delay / PATH;
    const double expected_half = ANGLE / 360.0 / 2.0 / speeds[i] / PATH;

    int found = 0;
    for (unsigned w = 0; w < windows.count; ++w) {
      const double centre = (windows.ranges[w].minimum + windows.ranges[w].maximum) / 2.0;
      if (fabs(centre - expected_centre) < 1e-9) {
        found = 1;
        CHECK_CLOSE((windows.ranges[w].maximum - windows.ranges[w].minimum) / 2.0,
                    expected_half, 1e-12);
      }
    }
    CHECK(found);
    if (windows.ranges) free(windows.ranges);
  }
}

/* A delay may be longer than one period; a phase could not.
 *
 * A phase wraps at a revolution, which quietly bounded the offset the library had to
 * cope with. A delay has no such bound -- a chopper can be set to open several periods
 * after the source -- and the rotation bracketing has to reach that far.
 */
static void test_a_delay_beyond_one_period_is_honoured(void) {
  TEST("a delay of several periods still opens a window");
  const double speed = 14.0;
  const double tau = 1.0 / speed;
  const double delay = 3.5 * tau;      /* three and a half revolutions after t=0 */

  const unsigned iv_bins = 1000;
  double * iv_edges = calloc(iv_bins + 1, sizeof(double));
  int * mask = calloc(iv_bins, sizeof(int));
  const chopper_parameters chopper = {speed, delay, ANGLE, PATH};
  const unsigned allowed = mask_one_chopper(chopper, mask, iv_edges, iv_bins);

  /* Openings recur at delay + n/speed. Within an arrival time of 0.1 s only n = -3
   * lands inside the grid, at 0.25 - 3/14 = 1/28 s, so at 1/280 s/m. */
  CHECK(allowed > 0);
  CHECK(allowed < iv_bins);

  unsigned counted = 0;
  const double centroid = allowed_centroid(mask, iv_bins, 1, iv_edges, &counted);
  CHECK_EQUAL_INT(counted, allowed);
  CHECK_CLOSE(centroid, (delay - 3.0 * tau) / PATH, 2.0 * 0.01 / (double) iv_bins);

  free(iv_edges);
  free(mask);
}

/* Which way the disk turns cannot change when a symmetric opening is on the path.
 *
 * The period is the same, the opening is centred on the delay either way, and a single
 * opening is symmetric about its own centre. So the two must give the same answer -- and
 * in particular the counter-rotating disk must still block something, rather than being
 * skipped and behaving as though it were permanently open.
 */
static void test_a_counter_rotating_disk_blocks_the_same_neutrons(void) {
  TEST("reversing a single-opening disk changes nothing");
  const unsigned iv_bins = 1000;
  double * forward_edges = calloc(iv_bins + 1, sizeof(double));
  double * reverse_edges = calloc(iv_bins + 1, sizeof(double));
  int * forward = calloc(iv_bins, sizeof(int));
  int * reverse = calloc(iv_bins, sizeof(int));

  const chopper_parameters cw  = {14.0, 0.02, ANGLE, PATH};
  const chopper_parameters ccw = {-14.0, 0.02, ANGLE, PATH};
  const unsigned forward_allowed = mask_one_chopper(cw, forward, forward_edges, iv_bins);
  const unsigned reverse_allowed = mask_one_chopper(ccw, reverse, reverse_edges, iv_bins);

  /* the chopper has to actually chop, or "identical" would be trivially true */
  CHECK(forward_allowed > 0);
  CHECK(forward_allowed < iv_bins);
  CHECK_EQUAL_INT(reverse_allowed, forward_allowed);

  int differing = 0;
  for (unsigned i = 0; i < iv_bins; ++i) if (forward[i] != reverse[i]) ++differing;
  CHECK_EQUAL_INT(differing, 0);

  free(forward_edges);
  free(reverse_edges);
  free(forward);
  free(reverse);
}

/* The openings a train admits are the ones every chopper admits. */
static void test_two_choppers_admit_only_their_overlap(void) {
  TEST("a chopper train admits only what all of its choppers admit");
  const chopper_parameters choppers[2] = {
    {14.0, 0.02, ANGLE, PATH},
    {14.0, 0.02, ANGLE, 2.0 * PATH},   /* twice as far, so half the inverse velocity */
  };

  double lower = 0.0, upper = 0.0;
  const unsigned windows = chopper_inverse_velocity_limits(
    &lower, &upper, 2, choppers, 0.0005, 0.01, 0.0);

  /* The first admits 0.02/10 = 0.002; the second admits 0.02/20 = 0.001 and, one period
   * later, 0.0914.../20 = 0.00457. Only the first chopper's n=1 opening at 0.00914 and
   * the second's are candidates -- nothing coincides, so the train passes nothing. */
  CHECK_EQUAL_INT(windows, 0);
}

/* The same two choppers, delayed so that openings do line up. */
static void test_a_train_passes_a_coincident_opening(void) {
  TEST("a chopper train passes the openings its choppers share");
  const double tau = 1.0 / 14.0;
  const chopper_parameters choppers[2] = {
    {14.0, 0.02, ANGLE, PATH},
    {14.0, 0.04, ANGLE, 2.0 * PATH},   /* twice as far, twice the delay: same 1/v */
  };

  range_set windows = chopper_inverse_velocity_windows(2, choppers, 0.0005, 0.01, 0.0);

  /* Both choppers see the same neutron at 0.002 s/m -- 10 m in 0.02 s, 20 m in 0.04 s --
   * and again one period later for the near chopper, two for the far one. The shared
   * opening is only as wide as the narrower of the two, which is the far chopper's:
   * the same angle covers half the inverse velocity across twice the flight path. */
  const double narrower_half = ANGLE / 360.0 / 2.0 / 14.0 / (2.0 * PATH);
  CHECK_EQUAL_INT(windows.count, 2);
  if (windows.count == 2) {
    CHECK_CLOSE((windows.ranges[0].minimum + windows.ranges[0].maximum) / 2.0,
                0.02 / PATH, 1e-12);
    CHECK_CLOSE((windows.ranges[1].minimum + windows.ranges[1].maximum) / 2.0,
                (0.02 + tau) / PATH, 1e-12);
    for (unsigned w = 0; w < 2; ++w) {
      CHECK_CLOSE((windows.ranges[w].maximum - windows.ranges[w].minimum) / 2.0,
                  narrower_half, 1e-12);
    }
  }
  if (windows.ranges) free(windows.ranges);
}

/* The reported limits envelope every window, including the gaps between them.
 *
 * `chopper_inverse_velocity_limits` returns a count precisely so a caller can tell that
 * its lower and upper bounds span inverse velocities the train does not actually pass.
 */
static void test_the_reported_limits_envelope_the_windows(void) {
  TEST("reported limits envelope every window a train passes");
  const double tau = 1.0 / 14.0;
  const chopper_parameters choppers[2] = {
    {14.0, 0.02, ANGLE, PATH},
    {14.0, 0.04, ANGLE, 2.0 * PATH},
  };
  const double narrower_half = ANGLE / 360.0 / 2.0 / 14.0 / (2.0 * PATH);

  double lower = 0.0, upper = 0.0;
  const unsigned windows = chopper_inverse_velocity_limits(
    &lower, &upper, 2, choppers, 0.0005, 0.01, 0.0);

  CHECK_EQUAL_INT(windows, 2);
  CHECK_CLOSE(lower, 0.02 / PATH - narrower_half, 1e-12);
  CHECK_CLOSE(upper, (0.02 + tau) / PATH + narrower_half, 1e-12);
}

int main(void) {
  test_the_delay_names_the_same_window_at_any_speed();
  test_a_delay_beyond_one_period_is_honoured();
  test_a_counter_rotating_disk_blocks_the_same_neutrons();
  test_two_choppers_admit_only_their_overlap();
  test_a_train_passes_a_coincident_opening();
  test_the_reported_limits_envelope_the_windows();
  return chopper_test_report();
}
