/* Disk choppers described the way the NeXus standard and the McStas component describe
 * them: a flat, increasing list of angles, two per opening.
 *
 * `chopper_inverse_velocity_windows` and its wavelength, envelope and mask wrappers all
 * read the same `edges` array, and all place an edge at angle `a` on the beam path at
 *
 *     t(a) = delay + (beam - a) / (360 * speed)
 *
 * Four things have to hold for that to describe a real disk, and the tests below are
 * grouped by which one they pin:
 *
 *   - `beam` is the angle that is on the beam path at `delay`, and only the difference
 *     `beam - a` matters, so moving the mark and moving the openings are the same move;
 *   - every opening on the disk has to contribute, at every rotation;
 *   - an opening has to be placed by its angle, which means by the *signed* speed, since
 *     a disk turning the other way reaches the same angle at a different time;
 *   - the openings are consecutive pairs of the array, in the order given.
 *
 * The third is invisible on an opening symmetric about the mark, which reversal maps
 * onto itself. It only shows up on a disk whose openings are asymmetric, which is
 * exactly the disk this description exists to allow.
 */
#include <stdlib.h>
#include <string.h>
#include "chopper-lib.h"
#include "test_util.h"

#if !defined(CHOPPER_LIB_VERSION) || CHOPPER_LIB_VERSION < 40000
#error "These tests describe disks by a flat edge array, which needs chopper-lib 4.0.0 or newer"
#endif

static const double PATH = 10.0;     /* m from the source to the chopper */
static const double SPEED = 14.0;    /* Hz */
static const double DELAY = 0.02;    /* s, when the angle `beam` is on the beam path */
static const double IV_MIN = 0.0005; /* s/m */
static const double IV_MAX = 0.01;   /* s/m */

/* The build passes these to chopper-lib.c as CHOPPER_LIB_DEFINITIONS, and McStas
 * supplies its own; a test of a unit conversion has to name the units it expects, so
 * repeat them here rather than take whatever the library was compiled with. */
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
 * Where the mark is, and where the openings are relative to it.
 * ------------------------------------------------------------------------------- */

/* An opening symmetric about the mark is centred on the delay, however the disk turns.
 *
 * This is the description a fixed-width chopper reduces to -- `{-width/2, +width/2}` --
 * and the one case where the direction of rotation cannot matter, because reflecting a
 * symmetric opening maps it onto itself. The four disks here are the awkward ones: a
 * reversed disk, a delay of several periods, and a fast, narrow, distant chopper.
 */
static void test_a_symmetric_opening_is_centred_on_the_delay(void) {
  TEST("an opening symmetric about the mark is centred on the delay");
  const double widths[4] = {3.6, 3.6, 12.0, 0.7};
  double edges[4][2];
  for (int i = 0; i < 4; ++i) { edges[i][0] = -widths[i] / 2.0; edges[i][1] = widths[i] / 2.0; }
  const chopper_parameters disks[4] = {
    { SPEED, DELAY,      0.0, 2, edges[0], PATH},
    {-SPEED, DELAY,      0.0, 2, edges[1], PATH},  /* reversed: symmetric, same answer */
    { SPEED, 3.5/SPEED,  0.0, 2, edges[2], PATH},  /* a delay of several periods */
    { 71.0, -0.004,      0.0, 2, edges[3], 162.0}, /* fast, narrow and far away */
  };

  for (int i = 0; i < 4; ++i) {
    range_set got = chopper_inverse_velocity_windows(1, &disks[i], IV_MIN, IV_MAX, 0.0);
    const double period = 1.0 / fabs(disks[i].speed) / disks[i].path;  /* in inverse velocity */
    const double half = window_half_width(widths[i], disks[i].speed, disks[i].path);
    const double reference = disks[i].delay / disks[i].path;

    CHECK(got.count > 0);   /* an empty answer would satisfy the loop trivially */
    for (unsigned w = 0; w < got.count; ++w) {
      CHECK_CLOSE(half_of(got.ranges[w]), half, 1e-15);
      /* every centre is a whole number of turns away from the delay */
      const double turns = (centre_of(got.ranges[w]) - reference) / period;
      CHECK_CLOSE(turns - floor(turns + 0.5), 0.0, 1e-9);
      /* and consecutive windows are exactly one turn apart */
      if (w) CHECK_CLOSE(centre_of(got.ranges[w]) - centre_of(got.ranges[w-1]), period, 1e-15);
    }
    if (got.ranges) free(got.ranges);
  }
}

/* Only `beam - a` matters, so turning the mark and turning the disk are the same turn.
 *
 * `beam` exists so a caller can hand over the component's own numbers: the openings in
 * the disk's own frame, and separately where the beam crosses it. Nothing downstream may
 * depend on the two independently, or that hand-over would not be faithful.
 */
static void test_the_beam_angle_and_the_openings_move_together(void) {
  TEST("moving the mark by an angle is moving every opening by it");
  double at_zero[2] = {0.0, 10.0};
  double moved[2] = {30.0, 40.0};
  const chopper_parameters disks[3] = {
    {SPEED, DELAY,  0.0, 2, at_zero, PATH},  /* openings at 0 and the mark at 0 */
    {SPEED, DELAY, 30.0, 2, moved,   PATH},  /* both moved by 30 degrees */
    {SPEED, DELAY, 30.0, 2, at_zero, PATH},  /* only the mark moved: a different disk */
  };

  range_set a = chopper_inverse_velocity_windows(1, &disks[0], IV_MIN, IV_MAX, 0.0);
  range_set b = chopper_inverse_velocity_windows(1, &disks[1], IV_MIN, IV_MAX, 0.0);
  range_set c = chopper_inverse_velocity_windows(1, &disks[2], IV_MIN, IV_MAX, 0.0);

  CHECK(a.count > 0);
  CHECK_EQUAL_INT(b.count, a.count);
  if (b.count == a.count) {
    for (unsigned w = 0; w < a.count; ++w) {
      CHECK_CLOSE(b.ranges[w].minimum, a.ranges[w].minimum, 1e-15);
      CHECK_CLOSE(b.ranges[w].maximum, a.ranges[w].maximum, 1e-15);
    }
  }
  /* moving the mark alone brings the openings onto the beam 30 degrees of a turn later */
  const double shift = 30.0 / 360.0 / SPEED / PATH;
  CHECK(c.count > 0);
  if (a.count && c.count) CHECK_CLOSE(centre_of(c.ranges[0]) - centre_of(a.ranges[0]), shift, 1e-15);

  if (a.ranges) free(a.ranges);
  if (b.ranges) free(b.ranges);
  if (c.ranges) free(c.ranges);
}

/* ---------------------------------------------------------------------------------
 * Every opening contributes, at every rotation.
 * ------------------------------------------------------------------------------- */

/* Six openings 60 degrees apart admit six times as many neutrons as one.
 *
 * The openings are evenly spaced, so they recur at a sixth of the disk period and the
 * answer is a comb six times as fine as the one-opening answer.
 *
 * This is also the case that used to run off the end of its own allocation: the range
 * array was sized for one range per rotation, but the loop fills one per rotation *per
 * opening*, so any disk with more than one opening wrote past the end of the buffer.
 */
static void test_every_opening_on_an_evenly_spaced_disk_admits_a_window(void) {
  TEST("each opening of an evenly spaced disk admits its own window");
  const unsigned n_openings = 6;
  const double width = 4.0;               /* degrees, per opening */
  double edges[12];
  for (unsigned i = 0; i < n_openings; ++i) {
    const double centre = 360.0 * (double) i / (double) n_openings;
    edges[2 * i] = centre - width / 2.0;
    edges[2 * i + 1] = centre + width / 2.0;
  }
  const chopper_parameters disk = {
    .speed = SPEED, .delay = DELAY, .beam = 0.0, .edge_count = 12, .edges = edges, .path = PATH
  };

  range_set got = chopper_inverse_velocity_windows(1, &disk, IV_MIN, IV_MAX, 0.0);

  /* Openings are on the beam at DELAY - k / (6 * SPEED) for every integer k, which over
   * all k is the same comb as DELAY + k / (6 * SPEED); in inverse velocity that is
   * (DELAY + k / 84) / PATH. Of those, k = -1 through 6 fall inside [IV_MIN, IV_MAX]
   * whole; k = -2 and k = 7 fall outside it entirely. */
  const double spacing = 1.0 / ((double) n_openings * SPEED) / PATH;
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
 * speed", so it cannot catch a fix that only multiplies the opening count. Two openings
 * a quarter turn apart can: the gaps they leave alternate between a quarter period and
 * three quarters of one.
 */
static void test_openings_are_placed_by_angle(void) {
  TEST("openings land where their angles put them, not evenly spaced");
  double edges[4] = {-1.0, 1.0, 89.0, 91.0};
  const chopper_parameters disk = {
    .speed = SPEED, .delay = DELAY, .beam = 0.0, .edge_count = 4, .edges = edges, .path = PATH
  };

  range_set got = chopper_inverse_velocity_windows(1, &disk, IV_MIN, IV_MAX, 0.0);

  /* A larger angle is reached *earlier*, so beam crossings are at DELAY + n / SPEED and
   * DELAY - (90/360) / SPEED + n / SPEED. Inside an arrival time of PATH * IV_MAX = 0.1 s
   * that is 0.02, 0.02 + 3/(4*14) and 0.02 + 1/14 -- the quarter-turn gap now falls after
   * the three-quarter one rather than before it. */
  const double expected[3] = {
    DELAY / PATH,
    (DELAY + 0.75 / SPEED) / PATH,
    (DELAY + 1.0 / SPEED) / PATH,
  };
  CHECK_EQUAL_INT(got.count, 3);
  if (got.count == 3) {
    for (unsigned w = 0; w < 3; ++w) CHECK_CLOSE(centre_of(got.ranges[w]), expected[w], 1e-15);
    /* the gaps really are uneven: three quarters of a turn, then a quarter */
    const double first_gap = centre_of(got.ranges[1]) - centre_of(got.ranges[0]);
    const double second_gap = centre_of(got.ranges[2]) - centre_of(got.ranges[1]);
    CHECK_CLOSE(first_gap / second_gap, 3.0, 1e-9);
  }
  if (got.ranges) free(got.ranges);
}

/* An opening is a consecutive pair of the array, and the array is in increasing order.
 *
 * Earlier versions took a list of {min, max} structures and were indifferent to the order
 * the openings arrived in, or to a pair given maximum-first. A flat array cannot be: the
 * pairing is positional, and the disk's angular extent is read off the first and last
 * entries. That is the cost of describing a disk the way NeXus and the component do, and
 * this test is here so the requirement is stated somewhere that fails when it stops
 * holding -- not to bless the answer that unsorted input happens to produce.
 */
static void test_openings_are_consecutive_pairs_in_increasing_order(void) {
  TEST("openings are consecutive pairs of the edge array, in the order given");
  double ordered[6] = {-2.0, 2.0, 58.0, 62.0, 178.0, 182.0};
  /* the same six angles, re-grouped: a caller who shuffles them describes another disk */
  double regrouped[6] = {-2.0, 58.0, 62.0, 178.0, 182.0, 2.0};
  const chopper_parameters a = {
    .speed = SPEED, .delay = DELAY, .beam = 0.0, .edge_count = 6, .edges = ordered, .path = PATH};
  const chopper_parameters b = {
    .speed = SPEED, .delay = DELAY, .beam = 0.0, .edge_count = 6, .edges = regrouped, .path = PATH};

  range_set want = chopper_inverse_velocity_windows(1, &a, IV_MIN, IV_MAX, 0.0);
  range_set got = chopper_inverse_velocity_windows(1, &b, IV_MIN, IV_MAX, 0.0);

  /* three openings, each 4 degrees wide, at 0, 60 and 180 degrees from the mark */
  CHECK(want.count > 0);
  const double half = window_half_width(4.0, SPEED, PATH);
  for (unsigned w = 0; w < want.count; ++w) CHECK_CLOSE(half_of(want.ranges[w]), half, 1e-15);

  /* the regrouped array is three 56, 116 and 180 degree openings: wider, and elsewhere */
  CHECK(got.count > 0);
  if (got.count) CHECK(half_of(got.ranges[0]) > 2.0 * half);

  if (want.ranges) free(want.ranges);
  if (got.ranges) free(got.ranges);
}

/* ---------------------------------------------------------------------------------
 * Which way the disk turns.
 * ------------------------------------------------------------------------------- */

/* Reversing the disk reflects its openings about the delay.
 *
 * `delay` fixes when the angle `beam` is on the beam path, whichever way the disk turns.
 * An opening ahead of that point in the direction of travel is reached before it; reverse
 * the disk and the same opening is reached after it. An opening centred on the mark
 * cannot show this -- reflecting it maps it onto itself -- so it takes an asymmetric one.
 */
static void test_reversing_the_disk_reflects_asymmetric_openings(void) {
  TEST("reversing the disk reflects its openings about the delay");
  double edges[2] = {0.0, 10.0};  /* entirely on one side of the mark */
  const chopper_parameters forward = {
    .speed = SPEED, .delay = DELAY, .beam = 0.0, .edge_count = 2, .edges = edges, .path = PATH};
  const chopper_parameters reverse = {
    .speed = -SPEED, .delay = DELAY, .beam = 0.0, .edge_count = 2, .edges = edges, .path = PATH};

  range_set f = chopper_inverse_velocity_windows(1, &forward, IV_MIN, IV_MAX, 0.0);
  range_set r = chopper_inverse_velocity_windows(1, &reverse, IV_MIN, IV_MAX, 0.0);

  /* forward: the 10 degree edge reaches the beam first and the mark arrives 10 degrees
   * later, so the opening runs up to the delay rather than away from it */
  const double reference = DELAY / PATH;
  const double span = 10.0 / 360.0 / SPEED / PATH;
  CHECK(f.count > 0);
  CHECK_EQUAL_INT(r.count, f.count);
  if (f.count > 0 && r.count == f.count) {
    CHECK_CLOSE(f.ranges[0].minimum, reference - span, 1e-15);
    CHECK_CLOSE(f.ranges[0].maximum, reference, 1e-15);
    /* reversed, the same opening follows the mark instead of preceding it */
    CHECK_CLOSE(r.ranges[0].minimum, reference, 1e-15);
    CHECK_CLOSE(r.ranges[0].maximum, reference + span, 1e-15);
  }
  if (f.ranges) free(f.ranges);
  if (r.ranges) free(r.ranges);
}

/* Where the mask puts its allowed bins, as contiguous runs of inverse velocity. */
static unsigned mask_runs(const chopper_parameters * disk, range * runs, unsigned max_runs) {
  const unsigned bins = 200000;
  double * grid = calloc(bins + 1, sizeof(double));
  int * mask = calloc(bins, sizeof(int));
  for (unsigned i = 0; i <= bins; ++i) grid[i] = IV_MAX * (double) i / (double) bins;
  const double times[2] = {0.0, 1.0e-12};   /* one vanishingly short emission bin */

  chopper_inverse_velocity_time_mask(
    mask, bins, 1, grid, bins + 1, times, 2, disk, 1, 0 /* no growing */);

  unsigned found = 0;
  double start = 0.0;
  int inside = 0;
  for (unsigned i = 0; i < bins; ++i) {
    if (mask[i] == CHOPPER_MASK_INCLUDED && !inside) { inside = 1; start = grid[i]; }
    else if (mask[i] != CHOPPER_MASK_INCLUDED && inside) {
      inside = 0;
      if (found < max_runs) { runs[found].minimum = start; runs[found].maximum = grid[i]; }
      ++found;
    }
  }
  if (inside && found < max_runs) { runs[found].minimum = start; runs[found].maximum = IV_MAX; ++found; }
  free(grid);
  free(mask);
  return found;
}

/* The mask and the window list describe the same chopper.
 *
 * They are separate implementations -- one intersects ranges, the other rejects bins --
 * so they can drift apart. They did: the mask placed opening angles with the unsigned
 * speed, which agrees for an opening symmetric about the mark and disagrees for any
 * other one, in the direction of mirroring the whole disk.
 */
static void test_the_mask_agrees_with_the_window_list(void) {
  TEST("the time mask admits what the window list says it should");
  double edges[4] = {0.0, 10.0, 100.0, 104.0};
  const double speeds[2] = {SPEED, -SPEED};

  for (int s = 0; s < 2; ++s) {
    const chopper_parameters disk = {
      .speed = speeds[s], .delay = DELAY, .beam = 0.0,
      .edge_count = 4, .edges = edges, .path = PATH};

    range_set want = chopper_inverse_velocity_windows(1, &disk, 0.0, IV_MAX, 0.0);
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

/* The mask moves with the delay, and repeats after exactly one period.
 *
 * `delay` says when the angle `beam` is on the beam path, so every window the mask
 * admits is fixed to it, and the openings recur every `1/|speed|` -- a delay a whole
 * period on describes the same disk and has to give the same mask back.
 *
 * Both halves are worth pinning here rather than leaving to the window list, which the
 * mask is only cross-checked against at one delay. Folding the delay into
 * `chopper_edge_time` left each of the two public functions binding `t0 = delay` and
 * never reading it again, which reads exactly like a mask that has stopped placing its
 * openings in time at all. Nothing in the suite said otherwise.
 */
static void test_the_mask_follows_the_delay(void) {
  TEST("the time mask moves with the delay and repeats after one period");
  double edges[2] = {-5.0, 5.0};      /* one opening, symmetric about the mark */
  const double tau = 1.0 / SPEED;
  /* the mask is generous by up to a bin at each end of a run, as above */
  const double tolerance = 3.0 * IV_MAX / 200000.0;

  /* A quarter period at a time across one whole turn. The opening is on the mark, so
   * it is on the beam path at `delay` itself, and every period from there: the mask
   * has to admit `delay / path`, and nothing that is not a whole number of turns from
   * it. A mask that ignored the delay would put its comb on the turns alone, which
   * DELAY -- 0.28 of a period -- is not on. */
  for (int q = 0; q < 4; ++q) {
    const double delay = DELAY + (double) q * tau / 4.0;
    const chopper_parameters disk = {
      .speed = SPEED, .delay = delay, .beam = 0.0,
      .edge_count = 2, .edges = edges, .path = PATH};
    range runs[8];
    const unsigned found = mask_runs(&disk, runs, 8);

    CHECK(found > 0);
    int admits_the_delay = 0;
    for (unsigned w = 0; w < found && w < 8; ++w) {
      const double centre = centre_of(runs[w]);
      if (fabs(centre - delay / PATH) <= tolerance) admits_the_delay = 1;
      /* and every window it does admit is a whole number of turns from that one */
      const double turns = (centre * PATH - delay) / tau;
      CHECK_CLOSE(turns, round(turns), tolerance * PATH / tau);
    }
    CHECK(admits_the_delay);
  }

  /* One period on is the same disk, bin for bin. */
  const chopper_parameters early = {
    .speed = SPEED, .delay = DELAY, .beam = 0.0,
    .edge_count = 2, .edges = edges, .path = PATH};
  const chopper_parameters late = {
    .speed = SPEED, .delay = DELAY + tau, .beam = 0.0,
    .edge_count = 2, .edges = edges, .path = PATH};
  range early_runs[8], late_runs[8];
  const unsigned early_found = mask_runs(&early, early_runs, 8);
  const unsigned late_found = mask_runs(&late, late_runs, 8);

  CHECK(early_found > 0);
  CHECK_EQUAL_INT(late_found, early_found);
  if (late_found == early_found) {
    for (unsigned w = 0; w < early_found && w < 8; ++w) {
      CHECK_CLOSE(late_runs[w].minimum, early_runs[w].minimum, tolerance);
      CHECK_CLOSE(late_runs[w].maximum, early_runs[w].maximum, tolerance);
    }
  }
}

/* ---------------------------------------------------------------------------------
 * Trains, envelopes and wavelengths.
 * ------------------------------------------------------------------------------- */

/* A parked disk is open or shut, and which one is a question about its angles.
 *
 * With no speed there is no period to recur on and no delay to apply, so all that is
 * left is whether the beam crosses an opening. Angles fold, which is what lets an
 * opening be written across the mark and a beam angle be given negative.
 */
static void test_a_parked_disk_is_open_only_when_the_beam_is_in_an_opening(void) {
  TEST("a parked disk stands open only when the beam crosses one of its openings");
  double edges[4] = {0.0, 10.0, 100.0, 104.0};
  double across_the_mark[2] = {350.0, 370.0};
  double about_the_mark[2] = {-85.0, 85.0};

  const chopper_parameters in_the_first = {
    .speed = 0.0, .delay = DELAY, .beam = 5.0, .edge_count = 4, .edges = edges, .path = PATH};
  const chopper_parameters in_the_second = {
    .speed = 0.0, .delay = DELAY, .beam = 102.0, .edge_count = 4, .edges = edges, .path = PATH};
  const chopper_parameters on_the_body = {
    .speed = 0.0, .delay = DELAY, .beam = 50.0, .edge_count = 4, .edges = edges, .path = PATH};
  /* on an edge: the opening runs from its first edge, and stops short of its last */
  const chopper_parameters on_the_opening_edge = {
    .speed = 0.0, .delay = DELAY, .beam = 0.0, .edge_count = 4, .edges = edges, .path = PATH};
  const chopper_parameters on_the_closing_edge = {
    .speed = 0.0, .delay = DELAY, .beam = 10.0, .edge_count = 4, .edges = edges, .path = PATH};
  /* an opening written across the mark, reached from either side of it */
  const chopper_parameters past_the_mark = {
    .speed = 0.0, .delay = DELAY, .beam = 5.0, .edge_count = 2, .edges = across_the_mark, .path = PATH};
  const chopper_parameters before_the_mark = {
    .speed = 0.0, .delay = DELAY, .beam = 355.0, .edge_count = 2, .edges = across_the_mark, .path = PATH};
  const chopper_parameters opposite_the_mark = {
    .speed = 0.0, .delay = DELAY, .beam = 180.0, .edge_count = 2, .edges = across_the_mark, .path = PATH};
  /* and one written about the mark with negative angles, which this library allows */
  const chopper_parameters inside_a_negative_opening = {
    .speed = 0.0, .delay = DELAY, .beam = -80.0, .edge_count = 2, .edges = about_the_mark, .path = PATH};
  const chopper_parameters outside_a_negative_opening = {
    .speed = 0.0, .delay = DELAY, .beam = 100.0, .edge_count = 2, .edges = about_the_mark, .path = PATH};
  /* a disk with no openings is solid, so it is shut wherever the beam crosses it */
  const chopper_parameters no_openings = {
    .speed = 0.0, .delay = DELAY, .beam = 0.0, .edge_count = 0, .edges = NULL, .path = PATH};

  CHECK_EQUAL_INT(chopper_parked_is_open(in_the_first), 1);
  CHECK_EQUAL_INT(chopper_parked_is_open(in_the_second), 1);
  CHECK_EQUAL_INT(chopper_parked_is_open(on_the_body), 0);
  CHECK_EQUAL_INT(chopper_parked_is_open(on_the_opening_edge), 1);
  CHECK_EQUAL_INT(chopper_parked_is_open(on_the_closing_edge), 0);
  CHECK_EQUAL_INT(chopper_parked_is_open(past_the_mark), 1);
  CHECK_EQUAL_INT(chopper_parked_is_open(before_the_mark), 1);
  CHECK_EQUAL_INT(chopper_parked_is_open(opposite_the_mark), 0);
  CHECK_EQUAL_INT(chopper_parked_is_open(inside_a_negative_opening), 1);
  CHECK_EQUAL_INT(chopper_parked_is_open(outside_a_negative_opening), 0);
  CHECK_EQUAL_INT(chopper_parked_is_open(no_openings), 0);
}

/* A disk parked shut empties the answer, and says which disk did it.
 *
 * It is a beam stop: nothing gets past it at any time, so the train admits no inverse
 * velocity and its mask has no open bin. That is the same answer a disk with no openings
 * gives, and the same one a train whose disks never agree gives, which is what makes it
 * the consistent one -- a caller already has to handle an empty result from either.
 *
 * What an empty result cannot say is *which* disk emptied it, and a disk parked shut is
 * nearly always one left mis-set. Hence the message, which names it. Watch for it in the
 * log; the checks here are on the answer, which is the part with a value to compare.
 */
static void test_a_disk_parked_shut_blocks_the_whole_train(void) {
  TEST("a disk parked with the beam on its body passes nothing, and is named");
  double edges[2] = {-2.0, 2.0};
  const chopper_parameters train[2] = {
    {.speed = SPEED, .delay = DELAY, .beam = 0.0, .edge_count = 2, .edges = edges, .path = PATH},
    /* parked with the beam a quarter turn from its only opening */
    {.speed = 0.0, .delay = DELAY, .beam = 90.0, .edge_count = 2, .edges = edges, .path = PATH},
  };
  CHECK_EQUAL_INT(chopper_parked_is_open(train[1]), 0);

  /* the turning disk on its own admits something, so the shut one is what empties it */
  range_set turning = chopper_inverse_velocity_windows(1, train, IV_MIN, IV_MAX, 0.0);
  range_set with_shut = chopper_inverse_velocity_windows(2, train, IV_MIN, IV_MAX, 0.0);
  CHECK(turning.count > 0);
  CHECK_EQUAL_INT(with_shut.count, 0);
  if (turning.ranges) free(turning.ranges);
  if (with_shut.ranges) free(with_shut.ranges);

  /* the envelope reports nothing and leaves the bounds alone, as it does for any train
   * that passes nothing -- the contract a caller already has to handle */
  double lower = -1.0, upper = -1.0;
  const unsigned count = chopper_inverse_velocity_limits(
    &lower, &upper, 2, train, IV_MIN, IV_MAX, 0.0);
  CHECK_EQUAL_INT(count, 0);
  CHECK_CLOSE(lower, -1.0, 0.0);
  CHECK_CLOSE(upper, -1.0, 0.0);

  /* and the mask masks off every bin, which is what the beam stop path does */
  const unsigned bins = 2000;
  double * grid = calloc(bins + 1, sizeof(double));
  int * alone = calloc(bins, sizeof(int));
  int * with_shut_mask = calloc(bins, sizeof(int));
  for (unsigned i = 0; i <= bins; ++i) grid[i] = IV_MAX * (double) i / (double) bins;
  const double times[2] = {0.0, 1.0e-12};

  const unsigned open_alone = chopper_inverse_velocity_time_mask(
    alone, bins, 1, grid, bins + 1, times, 2, train, 1, 0);
  const unsigned open_with_shut = chopper_inverse_velocity_time_mask(
    with_shut_mask, bins, 1, grid, bins + 1, times, 2, train, 2, 0);

  CHECK(open_alone > 0);
  CHECK_EQUAL_INT(open_with_shut, 0);
  int every_bin_masked = 1;
  for (unsigned i = 0; i < bins; ++i) if (with_shut_mask[i] == CHOPPER_MASK_INCLUDED) every_bin_masked = 0;
  CHECK(every_bin_masked);

  free(grid);
  free(alone);
  free(with_shut_mask);
}

/* Two openings offer twice as many chances to pass; a second disk takes half of them back.
 *
 * The two-opening disk admits a window every half period. The single-opening disk behind
 * it admits one every whole period, at inverse velocities the first disk also admits, so
 * the train passes every other one of the first disk's windows and none of its own extra.
 */
static void test_a_train_keeps_only_the_openings_every_disk_admits(void) {
  TEST("a train of multi-opening disks passes only their common openings");
  double two[4] = {-2.0, 2.0, 178.0, 182.0};
  double one[2] = {-1.0, 1.0};   /* half as wide */
  const chopper_parameters train[2] = {
    {.speed = SPEED, .delay = DELAY, .beam = 0.0, .edge_count = 4, .edges = two, .path = PATH},
    {.speed = SPEED, .delay = DELAY, .beam = 0.0, .edge_count = 2, .edges = one, .path = PATH},
  };

  range_set alone = chopper_inverse_velocity_windows(1, &train[0], IV_MIN, IV_MAX, 0.0);
  range_set both = chopper_inverse_velocity_windows(2, train, IV_MIN, IV_MAX, 0.0);

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

/* A disk parked open is not a closed one: the window functions step over it.
 *
 * Its opening is on the beam and stays there, so it passes every inverse velocity and
 * constrains nothing. A disk parked *shut* is the opposite of that, and is tested below.
 */
static void test_a_disk_parked_open_is_ignored(void) {
  TEST("a disk parked open does not block anything");
  double edges[2] = {-2.0, 2.0};
  const chopper_parameters train[2] = {
    {.speed = SPEED, .delay = DELAY, .beam = 0.0, .edge_count = 2, .edges = edges, .path = PATH},
    {.speed = 0.0, .delay = DELAY, .beam = 0.0, .edge_count = 2, .edges = edges, .path = PATH},
  };

  range_set turning = chopper_inverse_velocity_windows(1, train, IV_MIN, IV_MAX, 0.0);
  range_set with_still = chopper_inverse_velocity_windows(2, train, IV_MIN, IV_MAX, 0.0);

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
 * A multi-opening disk makes this the common case rather than the exception, which is why
 * the count comes back alongside the bounds: two windows and a gap between them report
 * the same pair of numbers as one window covering the whole span.
 */
static void test_the_envelope_spans_the_outermost_windows(void) {
  TEST("the reported envelope spans every window, including the gaps");
  double edges[4] = {-2.0, 2.0, 178.0, 182.0};
  const chopper_parameters disk = {
    .speed = SPEED, .delay = DELAY, .beam = 0.0, .edge_count = 4, .edges = edges, .path = PATH};

  range_set expected = chopper_inverse_velocity_windows(1, &disk, IV_MIN, IV_MAX, 0.0);
  double lower = 0.0, upper = 0.0;
  const unsigned count = chopper_inverse_velocity_limits(
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
  double a[2] = {-0.5, 0.5};
  double b[2] = {-0.5, 0.5};
  /* the second disk is a quarter period out of step with the first, at the same distance */
  const chopper_parameters train[2] = {
    {.speed = SPEED, .delay = DELAY, .beam = 0.0, .edge_count = 2, .edges = a, .path = PATH},
    {.speed = SPEED, .delay = DELAY + 0.25 / SPEED, .beam = 0.0,
     .edge_count = 2, .edges = b, .path = PATH},
  };

  double lower = -1.0, upper = -1.0;
  const unsigned count = chopper_inverse_velocity_limits(
    &lower, &upper, 2, train, IV_MIN, IV_MAX, 0.0);

  CHECK_EQUAL_INT(count, 0);
  /* the documented contract: the bounds are only set when the count is non-zero */
  CHECK_CLOSE(lower, -1.0, 0.0);
  CHECK_CLOSE(upper, -1.0, 0.0);
}

/* The wavelength wrapper is the inverse-velocity answer in different units. */
static void test_wavelength_limits_are_the_inverse_velocity_limits_converted(void) {
  TEST("the wavelength envelope is the inverse velocity envelope, converted");
  double edges[2] = {-2.0, 2.0};
  const chopper_parameters disk = {
    .speed = SPEED, .delay = DELAY, .beam = 0.0, .edge_count = 2, .edges = edges, .path = PATH};
  const double lambda_min = 1.0, lambda_max = 20.0;   /* angstrom */

  double lo = 0.0, hi = 0.0;
  const unsigned count = chopper_wavelength_limits(
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
  /* and the answer is in angstrom, not seconds per metre */
  CHECK(lo > 1.0 && hi < 20.0);

  /* and the same chopper asked in inverse velocity gives the bounds those came from */
  double iv_lo = 0.0, iv_hi = 0.0;
  const unsigned iv_count = chopper_inverse_velocity_limits(
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
  double edges[4] = {-2.0, 2.0, 118.0, 122.0};
  const chopper_parameters disk = {
    .speed = SPEED, .delay = DELAY, .beam = 0.0, .edge_count = 4, .edges = edges, .path = PATH};
  const double emission = 2.86e-3;   /* s, roughly an ESS pulse */

  range_set sharp = chopper_inverse_velocity_windows(1, &disk, IV_MIN, IV_MAX, 0.0);
  range_set spread = chopper_inverse_velocity_windows(1, &disk, IV_MIN, IV_MAX, emission);

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

/* A disk with no openings in it is a beam stop.
 *
 * It is a degenerate description rather than a real chopper, but the two implementations
 * have to make the same of it or a caller that builds edges dynamically gets one answer
 * from the envelope and the opposite one from the mask.
 */
static void test_a_disk_with_no_openings_passes_nothing(void) {
  TEST("a disk with no openings blocks everything, both ways of asking");
  const chopper_parameters blocked = {
    .speed = SPEED, .delay = DELAY, .beam = 0.0, .edge_count = 0, .edges = NULL, .path = PATH};

  range_set windows = chopper_inverse_velocity_windows(1, &blocked, IV_MIN, IV_MAX, 0.0);
  CHECK_EQUAL_INT(windows.count, 0);
  if (windows.ranges) free(windows.ranges);

  range runs[4];
  CHECK_EQUAL_INT(mask_runs(&blocked, runs, 4), 0);
}

int main(void) {
  test_a_symmetric_opening_is_centred_on_the_delay();
  test_the_beam_angle_and_the_openings_move_together();
  test_every_opening_on_an_evenly_spaced_disk_admits_a_window();
  test_openings_are_placed_by_angle();
  test_openings_are_consecutive_pairs_in_increasing_order();
  test_reversing_the_disk_reflects_asymmetric_openings();
  test_the_mask_agrees_with_the_window_list();
  test_the_mask_follows_the_delay();
  test_a_train_keeps_only_the_openings_every_disk_admits();
  test_a_disk_parked_open_is_ignored();
  test_a_parked_disk_is_open_only_when_the_beam_is_in_an_opening();
  test_a_disk_parked_shut_blocks_the_whole_train();
  test_the_envelope_spans_the_outermost_windows();
  test_a_train_that_admits_nothing_reports_nothing();
  test_wavelength_limits_are_the_inverse_velocity_limits_converted();
  test_a_longer_pulse_widens_every_window_downwards();
  test_a_disk_with_no_openings_passes_nothing();
  return chopper_test_report();
}
