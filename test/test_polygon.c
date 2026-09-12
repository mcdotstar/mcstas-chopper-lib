/* The transmitted region as polygons.
 *
 * The claim these hold the code to is that the polygons are *exact*, so most of them
 * compare against something that shares no code with them: a point-in-set test worked out
 * from the disk definitions directly, the mask, and closed-form areas.
 */
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "chopper-lib.h"
#include "test_util.h"

/* ---- the BIFROST train, as chopcal sets it for a 3 AA band, beam apertures and all ---- */
static double ps1_edges[2] = {-85.0, 85.0};
static double ps2_edges[2] = {-85.0, 85.0};
static double fo1_edges[2] = {-19.13, 19.13};
static double fo2_edges[2] = {-26.005, 26.005};
static double bw1_edges[2] = {-80.5, 80.5};
static double bw2_edges[2] = {-80.5, 80.5};

static chopper_parameters bifrost[6] = {
  {196.0, 0.0038781004374078398, 0.0, 2, ps1_edges,  6.3420000000000005,  5.650725436158406},
  {196.0, 0.0060873974895620350, 0.0, 2, ps2_edges,  6.3620000000000000,  5.650725436158406},
  { 14.0, 0.0061399701448323090, 0.0, 2, fo1_edges,  8.5300000000000000,  5.650725436158406},
  { 14.0, 0.0095476374535256970, 0.0, 2, fo2_edges, 14.9730000000000000,  5.650725436158406},
  { 14.0, 0.0428822715471184300, 0.0, 2, bw1_edges, 78.0000000000000000, 13.379640178419610},
  {-14.0, 0.0428822715471184300, 0.0, 2, bw2_edges, 78.0200000000000000, 13.379640178419610},
};

#define SOURCE_INVERSE_VELOCITY_MIN 1.0e-4
#define SOURCE_INVERSE_VELOCITY_RANGE 1.13e-2   /* out to about 45 AA */
#define SOURCE_TIME_MIN 0.0
#define SOURCE_TIME_RANGE 3.0e-3

static chopper_polygon_set bifrost_source(void) {
  chopper_polygon_set set = chopper_polygon_set_empty();
  const chopper_polygon rectangle = chopper_polygon_rectangle(
      SOURCE_INVERSE_VELOCITY_MIN, SOURCE_INVERSE_VELOCITY_RANGE,
      SOURCE_TIME_MIN, SOURCE_TIME_RANGE);
  chopper_polygon_set_add(&set, &rectangle);
  return set;
}

/* ---- an independent answer: does this one neutron get through? -------------------- */
/* Worked from the disk definitions rather than from any polygon, so agreement between
 * the two is evidence rather than a tautology. */
static int neutron_passes(const chopper_parameters * choppers, const unsigned count,
                          const double inverse_velocity, const double time,
                          const double path_spread) {
  for (unsigned i = 0; i < count; ++i) {
    const chopper_parameters c = choppers[i];
    if (c.speed == 0.0) continue;
    const double tau = 1.0 / fabs(c.speed);
    const double aperture = c.aperture > 0 ? c.aperture / 2.0 / 360.0 / fabs(c.speed) : 0.0;
    const double early = time + c.path * inverse_velocity;
    const double late = time + (c.path + path_spread) * inverse_velocity;
    int through = 0;
    for (unsigned w = 0; w + 1 < c.edge_count && !through; w += 2) {
      double lo = c.delay + (c.beam - c.edges[w]) / 360.0 / c.speed;
      double hi = c.delay + (c.beam - c.edges[w + 1]) / 360.0 / c.speed;
      if (hi < lo) { const double s = lo; lo = hi; hi = s; }
      lo -= aperture; hi += aperture;
      const long first = (long) floor((early - hi) / tau) - 1;
      const long last = (long) ceil((late - lo) / tau) + 1;
      for (long n = first; n <= last && !through; ++n) {
        if (early <= hi + (double) n * tau && late >= lo + (double) n * tau) through = 1;
      }
    }
    if (!through) return 0;
  }
  return 1;
}

/* A deterministic low-discrepancy sequence, so the checks below do not depend on which
 * rand() the platform ships. */
static double radical_inverse(unsigned i, const unsigned base) {
  double f = 1.0 / (double) base, result = 0.0;
  while (i) { result += f * (double) (i % base); i /= base; f /= (double) base; }
  return result;
}

int main(void) {
  TEST("a rectangle is four vertices and the area it should be");
  {
    const chopper_polygon r = chopper_polygon_rectangle(2.0, 3.0, 10.0, 5.0);
    CHECK_EQUAL_INT(r.count, 4);
    CHECK_CLOSE(chopper_polygon_area(&r), 15.0, 1e-12);
    double lo = 0, hi = 0;
    chopper_polygon_extent(&r, 1.0, 0.0, &lo, &hi);
    CHECK_CLOSE(lo, 2.0, 1e-12);
    CHECK_CLOSE(hi, 5.0, 1e-12);
    chopper_polygon_extent(&r, 0.0, 1.0, &lo, &hi);
    CHECK_CLOSE(lo, 10.0, 1e-12);
    CHECK_CLOSE(hi, 15.0, 1e-12);
  }
  {
    const chopper_polygon none = chopper_polygon_rectangle(2.0, -1.0, 10.0, 5.0);
    CHECK_EQUAL_INT(none.count, 0);
    CHECK_CLOSE(chopper_polygon_area(&none), 0.0, 0.0);
  }

  TEST("clipping a half-plane keeps the half it is asked for");
  {
    chopper_polygon r = chopper_polygon_rectangle(0.0, 4.0, 0.0, 2.0);   /* area 8 */
    CHECK_EQUAL_INT(chopper_polygon_clip_halfplane(&r, 1.0, 0.0, 2.0), 1);
    CHECK_CLOSE(chopper_polygon_area(&r), 4.0, 1e-12);
    double hi = 0;
    chopper_polygon_extent(&r, 1.0, 0.0, NULL, &hi);
    CHECK_CLOSE(hi, 2.0, 1e-12);
  }
  {
    /* A diagonal cut through a unit square takes exactly half of it, and leaves a
     * triangle -- one vertex fewer, not more. */
    chopper_polygon r = chopper_polygon_rectangle(0.0, 1.0, 0.0, 1.0);
    CHECK_EQUAL_INT(chopper_polygon_clip_halfplane(&r, 1.0, 1.0, 1.0), 1);
    CHECK_CLOSE(chopper_polygon_area(&r), 0.5, 1e-12);
    CHECK_EQUAL_INT(r.count, 3);
  }
  {
    chopper_polygon r = chopper_polygon_rectangle(0.0, 1.0, 0.0, 1.0);
    CHECK_EQUAL_INT(chopper_polygon_clip_halfplane(&r, 1.0, 0.0, -1.0), 1);
    CHECK_EQUAL_INT(r.count, 0);       /* wholly outside */
    CHECK_CLOSE(chopper_polygon_area(&r), 0.0, 0.0);
  }
  {
    /* Clipping along an edge the polygon already has must not invent vertices. */
    chopper_polygon r = chopper_polygon_rectangle(0.0, 1.0, 0.0, 1.0);
    CHECK_EQUAL_INT(chopper_polygon_clip_halfplane(&r, 1.0, 0.0, 1.0), 1);
    CHECK_EQUAL_INT(r.count, 4);
    CHECK_CLOSE(chopper_polygon_area(&r), 1.0, 1e-12);
  }

  TEST("a wedge with equal paths is an ordinary slab");
  {
    chopper_polygon slab = chopper_polygon_rectangle(0.0, 1.0, 0.0, 1.0);
    chopper_polygon wedge = slab;
    CHECK_EQUAL_INT(chopper_polygon_clip_wedge(&slab, 0.5, 0.5, 1.0, 0.2, 0.8), 1);
    CHECK_EQUAL_INT(chopper_polygon_clip_wedge(&wedge, 0.5, 0.5, 1.0, 0.2, 0.8), 1);
    CHECK_CLOSE(chopper_polygon_area(&slab), chopper_polygon_area(&wedge), 0.0);
    /* Not 0.6: the lower line t = 0.2 - a/2 leaves the square at a = 0.4, so below that
     * the band is cut by the square's own floor.
     *   int_0^0.4 0.6 da  +  int_0.4^1 (0.8 - a/2) da  =  0.24 + 0.27 */
    CHECK_CLOSE(chopper_polygon_area(&slab), 0.51, 1e-12);
  }
  {
    /* Widening the far path can only add area, never remove it. */
    chopper_polygon narrow = chopper_polygon_rectangle(0.0, 1.0, 0.0, 1.0);
    chopper_polygon wide = narrow;
    chopper_polygon_clip_wedge(&narrow, 0.5, 0.5, 1.0, 0.2, 0.8);
    chopper_polygon_clip_wedge(&wide, 0.5, 0.7, 1.0, 0.2, 0.8);
    CHECK(chopper_polygon_area(&wide) > chopper_polygon_area(&narrow));
  }

  TEST("overflowing the vertex cap is refused, not truncated");
  {
    /* A regular polygon at the cap, cut so that exactly one vertex falls outside. That
     * keeps 31 and adds two crossings, which is one more than there is room for. Dropping
     * a vertex instead would quietly enlarge the polygon, so this must fail loudly. */
    chopper_polygon full;
    full.count = CHOPPER_POLYGON_MAX_VERTICES;
    double highest = -2.0;
    for (unsigned i = 0; i < full.count; ++i) {
      const double angle = 2.0 * 3.14159265358979323846 * (double) i / (double) full.count;
      full.vertex[i].inverse_velocity = cos(angle);
      full.vertex[i].time = sin(angle);
      if (full.vertex[i].time > highest) highest = full.vertex[i].time;
    }
    const chopper_polygon before = full;
    double second = -2.0;
    for (unsigned i = 0; i < full.count; ++i) {
      if (full.vertex[i].time < highest && full.vertex[i].time > second)
        second = full.vertex[i].time;
    }
    CHECK_EQUAL_INT(chopper_polygon_clip_halfplane(&full, 0.0, 1.0,
                                                   (highest + second) / 2.0), 0);
    CHECK_EQUAL_INT(full.count, before.count);
    CHECK_CLOSE(full.vertex[0].inverse_velocity, before.vertex[0].inverse_velocity, 0.0);
  }

  TEST("a disk that would overflow a polygon fails the whole transmission");
  {
    /* A source region need not be a rectangle -- the vertex budget is V + 2n, so a caller
     * may hand in a polygon that is already near the cap. One that reaches it, cut by a
     * disk that trims a single corner, cannot be represented and must not be truncated. */
    chopper_polygon ring;
    ring.count = CHOPPER_POLYGON_MAX_VERTICES;
    for (unsigned i = 0; i < ring.count; ++i) {
      const double angle = 2.0 * 3.14159265358979323846 * (double) i / (double) ring.count;
      ring.vertex[i].inverse_velocity = 1.0e-3 + 5.0e-4 * cos(angle);
      ring.vertex[i].time = 0.5 + 0.1 * sin(angle);       /* topmost vertex at 0.6 */
    }
    /* At the source, so the window is horizontal in time: 1 Hz and a 180 degree opening
     * is a half-second window inside a one-second turn, placed to cut just under 0.6 --
     * above the second-highest vertex at 0.5 + 0.1*cos(2*pi/32) = 0.598. */
    static double edges[2] = {-90.0, 90.0};
    const double half_window = 90.0 / 360.0 / 1.0;
    const chopper_parameters trimmer = {1.0, 0.599 - half_window, 0.0, 2, edges, 0.0, 0.0};

    chopper_polygon_set set = chopper_polygon_set_empty();
    CHECK_EQUAL_INT(chopper_polygon_set_add(&set, &ring), 1);
    CHECK_EQUAL_INT(set.count, 1);
    CHECK_EQUAL_INT(chopper_polygon_set_transmit(&set, trimmer, 0.0), 0);
    CHECK_EQUAL_INT(set.count, 0);        /* and it owns nothing after failing */
    CHECK(set.polygon == NULL);
    chopper_polygon_set_free(&set);
  }

  TEST("a negative path spread is refused");
  {
    chopper_polygon_set set = bifrost_source();
    CHECK_EQUAL_INT(chopper_polygon_set_transmit(&set, bifrost[0], -1.0), 0);
    CHECK_EQUAL_INT(set.count, 0);
    chopper_polygon_set_free(&set);
  }

  TEST("a disk with no openings stops everything");
  {
    chopper_parameters blank = bifrost[0];
    blank.edge_count = 0;
    blank.edges = NULL;
    chopper_polygon_set set = bifrost_source();
    CHECK_EQUAL_INT(chopper_polygon_set_transmit(&set, blank, 0.0), 1);
    CHECK_EQUAL_INT(set.count, 0);
    chopper_polygon_set_free(&set);
  }

  TEST("a set owns one allocation and drops what carries nothing");
  {
    chopper_polygon_set set = chopper_polygon_set_empty();
    CHECK_EQUAL_INT(set.count, 0);
    CHECK(set.polygon == NULL);
    const chopper_polygon good = chopper_polygon_rectangle(0.0, 1.0, 0.0, 1.0);
    for (unsigned i = 0; i < 40; ++i) CHECK_EQUAL_INT(chopper_polygon_set_add(&set, &good), 1);
    CHECK_EQUAL_INT(set.count, 40);
    CHECK_CLOSE(chopper_polygon_set_area(&set), 40.0, 1e-12);
    /* a line, not a region */
    chopper_polygon flat = chopper_polygon_rectangle(0.0, 1.0, 0.0, 1.0);
    chopper_polygon_clip_halfplane(&flat, 1.0, 0.0, 0.0);
    CHECK_EQUAL_INT(chopper_polygon_set_add(&set, &flat), 1);
    CHECK_EQUAL_INT(set.count, 40);
    chopper_polygon_set_free(&set);
    CHECK_EQUAL_INT(set.count, 0);
    CHECK(set.polygon == NULL);
    chopper_polygon_set_free(&set);   /* twice is fine */
    chopper_polygon_set_free(NULL);
  }

  TEST("one disk passes what it should, against a closed-form area");
  {
    /* A disk of one opening, at the source, so there is no shear: the accepted region is
     * the emission window crossed with the times the opening is on the beam. The opening
     * is 36 degrees at 10 Hz, so it is open 0.01 s in every 0.1 s, centred on delay. */
    static double edges[2] = {-18.0, 18.0};
    const chopper_parameters disk = {10.0, 0.005, 0.0, 2, edges, 0.0, 0.0};
    chopper_polygon_set set = chopper_polygon_set_empty();
    const chopper_polygon source = chopper_polygon_rectangle(1.0e-4, 1.0e-3, 0.0, 0.1);
    chopper_polygon_set_add(&set, &source);
    CHECK_EQUAL_INT(chopper_polygon_set_transmit(&set, disk, 0.0), 1);
    /* one whole opening inside [0, 0.1] */
    CHECK_CLOSE(chopper_polygon_set_area(&set), 1.0e-3 * 0.01, 1e-15);
    chopper_polygon_set_free(&set);
  }

  TEST("a disk parked open constrains nothing, one parked shut stops everything");
  {
    static double edges[2] = {-10.0, 10.0};
    chopper_parameters parked = {0.0, 0.0, 0.0, 2, edges, 1.0, 0.0};
    chopper_polygon_set set = chopper_polygon_set_empty();
    const chopper_polygon source = chopper_polygon_rectangle(1.0e-4, 1.0e-3, 0.0, 3.0e-3);
    chopper_polygon_set_add(&set, &source);
    const double before = chopper_polygon_set_area(&set);
    CHECK_EQUAL_INT(chopper_polygon_set_transmit(&set, parked, 0.0), 1);
    CHECK_CLOSE(chopper_polygon_set_area(&set), before, 0.0);

    parked.beam = 180.0;                       /* now the beam is on the disk body */
    CHECK_EQUAL_INT(chopper_polygon_set_transmit(&set, parked, 0.0), 1);
    CHECK_EQUAL_INT(set.count, 0);
    chopper_polygon_set_free(&set);
  }

  TEST("a disk open for a whole turn constrains nothing");
  {
    static double edges[2] = {-180.0, 180.0};
    const chopper_parameters always = {10.0, 0.0, 0.0, 2, edges, 5.0, 0.0};
    chopper_polygon_set set = chopper_polygon_set_empty();
    const chopper_polygon source = chopper_polygon_rectangle(1.0e-4, 1.0e-3, 0.0, 3.0e-3);
    chopper_polygon_set_add(&set, &source);
    const double before = chopper_polygon_set_area(&set);
    CHECK_EQUAL_INT(chopper_polygon_set_transmit(&set, always, 0.0), 1);
    CHECK_CLOSE(chopper_polygon_set_area(&set), before, 0.0);
    chopper_polygon_set_free(&set);
  }

  TEST("the BIFROST train: every polygon vertex is a neutron that gets through");
  {
    chopper_polygon_set set = bifrost_source();
    CHECK_EQUAL_INT(chopper_polygon_set_transmit_train(&set, 6, bifrost, NULL), 1);
    CHECK(set.count > 0);
    CHECK(chopper_polygon_set_area(&set) > 0.0);
    for (unsigned i = 0; i < set.count; ++i) {
      CHECK(set.polygon[i].count >= 3);
      CHECK(set.polygon[i].count <= CHOPPER_POLYGON_MAX_VERTICES);
      /* tighter than the cap: a rectangle plus two clips per disk */
      CHECK(set.polygon[i].count <= 4 + 2 * 6);
      /* The centroid is strictly inside a convex polygon, so it must pass. A vertex sits
       * exactly on a boundary, where a rounding step either way decides it. */
      double a = 0.0, t = 0.0;
      for (unsigned v = 0; v < set.polygon[i].count; ++v) {
        a += set.polygon[i].vertex[v].inverse_velocity;
        t += set.polygon[i].vertex[v].time;
      }
      a /= (double) set.polygon[i].count;
      t /= (double) set.polygon[i].count;
      CHECK(neutron_passes(bifrost, 6, a, t, 0.0));
    }
    chopper_polygon_set_free(&set);
  }

  TEST("the transmitted area agrees with counting points that get through");
  {
    chopper_polygon_set set = bifrost_source();
    chopper_polygon_set_transmit_train(&set, 6, bifrost, NULL);
    const double source_area = SOURCE_INVERSE_VELOCITY_RANGE * SOURCE_TIME_RANGE;
    const double from_polygons = chopper_polygon_set_area(&set) / source_area;

    const unsigned samples = 200000;
    unsigned through = 0;
    for (unsigned i = 1; i <= samples; ++i) {
      const double a = SOURCE_INVERSE_VELOCITY_MIN
                     + SOURCE_INVERSE_VELOCITY_RANGE * radical_inverse(i, 2);
      const double t = SOURCE_TIME_MIN + SOURCE_TIME_RANGE * radical_inverse(i, 3);
      if (neutron_passes(bifrost, 6, a, t, 0.0)) ++through;
    }
    const double counted = (double) through / (double) samples;
    /* A quasi-random sequence converges faster than 1/sqrt(N), but the boundary is a
     * thin diagonal band, so keep the tolerance honest at a few percent of the answer. */
    CHECK_CLOSE(from_polygons, counted, 0.05 * counted);
    chopper_polygon_set_free(&set);
  }

  TEST("the polygons and the mask agree on the band, where the window function does not");
  {
    chopper_polygon_set set = bifrost_source();
    chopper_polygon_set_transmit_train(&set, 6, bifrost, NULL);
    range_set bands = chopper_polygon_set_inverse_velocity_ranges(&set);
    CHECK_EQUAL_INT(bands.count, 1);

    /* The window function reports that band and a second one the train does not pass.
     * That is what this whole section exists to fix, so record it rather than assert it
     * away: if it ever stops being true, this check fails and should be deleted. */
    const double search_low = SOURCE_INVERSE_VELOCITY_MIN;
    const double search_high = SOURCE_INVERSE_VELOCITY_MIN + SOURCE_INVERSE_VELOCITY_RANGE;
    range_set windows = chopper_inverse_velocity_windows(6, bifrost, search_low,
                                                         search_high, SOURCE_TIME_RANGE);
    CHECK(windows.count > bands.count);

    /* Every polygon band must lie inside some window: the window function over-reports,
     * so it can never be the narrower of the two. */
    if (bands.count && windows.count) {
      int contained = 0;
      for (unsigned w = 0; w < windows.count; ++w) {
        if (windows.ranges[w].minimum <= bands.ranges[0].minimum * (1 + 1e-9)
            && windows.ranges[w].maximum >= bands.ranges[0].maximum * (1 - 1e-9)) contained = 1;
      }
      CHECK(contained);
    }
    if (bands.ranges) free(bands.ranges);
    if (windows.ranges) free(windows.ranges);
    chopper_polygon_set_free(&set);
  }

  TEST("the answer does not depend on the units it is computed in");
  {
    /* Inverse velocity in ms/km and time in ms rather than s/m and s: the two axes scale
     * by 1e6 and 1e3, so an area scales by 1e9. The transmitted *fraction* may not move.
     * A tolerance measured in absolute area would fail this. */
    chopper_polygon_set si = bifrost_source();
    chopper_polygon_set_transmit_train(&si, 6, bifrost, NULL);
    const double si_fraction = chopper_polygon_set_area(&si)
                             / (SOURCE_INVERSE_VELOCITY_RANGE * SOURCE_TIME_RANGE);

    chopper_parameters rescaled[6];
    for (unsigned i = 0; i < 6; ++i) {
      rescaled[i] = bifrost[i];
      rescaled[i].speed = bifrost[i].speed / 1.0e3;    /* kHz, so a period is in ms */
      rescaled[i].delay = bifrost[i].delay * 1.0e3;    /* ms */
      rescaled[i].path = bifrost[i].path / 1.0e3;      /* km */
    }
    chopper_polygon_set alt = chopper_polygon_set_empty();
    const chopper_polygon rectangle = chopper_polygon_rectangle(
        SOURCE_INVERSE_VELOCITY_MIN * 1.0e6, SOURCE_INVERSE_VELOCITY_RANGE * 1.0e6,
        SOURCE_TIME_MIN * 1.0e3, SOURCE_TIME_RANGE * 1.0e3);
    chopper_polygon_set_add(&alt, &rectangle);
    chopper_polygon_set_transmit_train(&alt, 6, rescaled, NULL);
    const double alt_fraction = chopper_polygon_set_area(&alt)
                              / (SOURCE_INVERSE_VELOCITY_RANGE * 1.0e6
                                 * SOURCE_TIME_RANGE * 1.0e3);

    CHECK_EQUAL_INT(alt.count, si.count);
    CHECK_CLOSE(alt_fraction, si_fraction, 1e-12 * si_fraction);
    chopper_polygon_set_free(&si);
    chopper_polygon_set_free(&alt);
  }

  TEST("a path spread only ever widens the answer, and only in time");
  {
    double previous = 0.0;
    for (unsigned step = 0; step < 4; ++step) {
      const double fraction = (double) step * 1.0e-4;
      double spreads[6];
      for (unsigned i = 0; i < 6; ++i) spreads[i] = fraction * bifrost[i].path;
      chopper_polygon_set set = bifrost_source();
      CHECK_EQUAL_INT(chopper_polygon_set_transmit_train(&set, 6, bifrost, spreads), 1);
      const double area = chopper_polygon_set_area(&set);
      CHECK(area >= previous);
      previous = area;
      for (unsigned i = 0; i < set.count; ++i) {
        double a = 0.0, t = 0.0;
        for (unsigned v = 0; v < set.polygon[i].count; ++v) {
          a += set.polygon[i].vertex[v].inverse_velocity;
          t += set.polygon[i].vertex[v].time;
        }
        a /= (double) set.polygon[i].count;
        t /= (double) set.polygon[i].count;
        CHECK(neutron_passes(bifrost, 6, a, t, fraction * 0.0));
      }
      chopper_polygon_set_free(&set);
    }
  }

  TEST("a path spread wide enough to overlap turns is refused, not miscounted");
  {
    double spreads[6];
    for (unsigned i = 0; i < 6; ++i) spreads[i] = 0.5 * bifrost[i].path;  /* absurd */
    chopper_polygon_set set = bifrost_source();
    CHECK_EQUAL_INT(chopper_polygon_set_transmit_train(&set, 6, bifrost, spreads), 0);
    CHECK_EQUAL_INT(set.count, 0);
    chopper_polygon_set_free(&set);
  }

  TEST("the sampler draws inside the region, and its acceptance is the area ratio");
  {
    chopper_polygon_set set = bifrost_source();
    chopper_polygon_set_transmit_train(&set, 6, bifrost, NULL);
    const double source_area = SOURCE_INVERSE_VELOCITY_RANGE * SOURCE_TIME_RANGE;
    chopper_polygon_sampler sampler = chopper_polygon_sampler_make(&set, source_area);
    CHECK(sampler.count > 0);
    CHECK_CLOSE(sampler.acceptance, chopper_polygon_set_area(&set) / source_area, 1e-15);

    const unsigned draws = 20000;
    unsigned outside = 0;
    for (unsigned i = 1; i <= draws; ++i) {
      double a = 0, t = 0;
      chopper_polygon_sampler_draw(&sampler, radical_inverse(i, 2), radical_inverse(i, 3),
                                   radical_inverse(i, 5), &a, &t);
      if (!neutron_passes(bifrost, 6, a, t, 0.0)) ++outside;
    }
    CHECK_EQUAL_INT(outside, 0);

    /* Uniform, so an equal-area half of the region takes an equal share of the draws.
     * Split the plane at the midpoint in inverse velocity and compare the two. */
    double low = 0, high = 0;
    chopper_polygon_extent(&set.polygon[0], 1.0, 0.0, &low, &high);
    for (unsigned i = 1; i < set.count; ++i) {
      double l = 0, h = 0;
      chopper_polygon_extent(&set.polygon[i], 1.0, 0.0, &l, &h);
      if (l < low) low = l;
      if (h > high) high = h;
    }
    const double middle = (low + high) / 2.0;
    chopper_polygon_set lower_half = chopper_polygon_set_empty();
    for (unsigned i = 0; i < set.count; ++i) {
      chopper_polygon piece = set.polygon[i];
      chopper_polygon_clip_halfplane(&piece, 1.0, 0.0, middle);
      chopper_polygon_set_add(&lower_half, &piece);
    }
    const double expected = chopper_polygon_set_area(&lower_half)
                          / chopper_polygon_set_area(&set);
    unsigned below = 0;
    for (unsigned i = 1; i <= draws; ++i) {
      double a = 0, t = 0;
      chopper_polygon_sampler_draw(&sampler, radical_inverse(i, 2), radical_inverse(i, 3),
                                   radical_inverse(i, 5), &a, &t);
      if (a < middle) ++below;
    }
    CHECK_CLOSE((double) below / (double) draws, expected, 0.02);

    chopper_polygon_set_free(&lower_half);
    chopper_polygon_sampler_free(&sampler);
    chopper_polygon_set_free(&set);
  }

  TEST("containment: inside, outside, and exactly on the boundary");
  {
    const chopper_polygon square = chopper_polygon_rectangle(0.0, 2.0, 0.0, 2.0);
    CHECK_EQUAL_INT(chopper_polygon_contains(&square, 1.0, 1.0), 1);
    CHECK_EQUAL_INT(chopper_polygon_contains(&square, 3.0, 1.0), 0);
    CHECK_EQUAL_INT(chopper_polygon_contains(&square, 1.0, -0.5), 0);
    CHECK_EQUAL_INT(chopper_polygon_contains(&square, 0.0, 1.0), 1);   /* on an edge */
    CHECK_EQUAL_INT(chopper_polygon_contains(&square, 2.0, 2.0), 1);   /* on a corner */

    /* Winding must not matter: the same square listed the other way round. */
    chopper_polygon reversed = square;
    for (unsigned i = 0; i < square.count; ++i)
      reversed.vertex[i] = square.vertex[square.count - 1 - i];
    CHECK_EQUAL_INT(chopper_polygon_contains(&reversed, 1.0, 1.0), 1);
    CHECK_EQUAL_INT(chopper_polygon_contains(&reversed, 3.0, 1.0), 0);

    const chopper_polygon empty = chopper_polygon_rectangle(0.0, -1.0, 0.0, 1.0);
    CHECK_EQUAL_INT(chopper_polygon_contains(&empty, 0.0, 0.0), 0);
    CHECK_EQUAL_INT(chopper_polygon_contains(NULL, 0.0, 0.0), 0);
  }

  TEST("a set contains what its polygons do, and the region agrees with the disks");
  {
    chopper_polygon_set set = bifrost_source();
    chopper_polygon_set_transmit_train(&set, 6, bifrost, NULL);
    CHECK_EQUAL_INT(chopper_polygon_set_contains(NULL, 0.0, 0.0), 0);

    /* Centroids are inside, and the point test agrees a neutron there gets through. */
    for (unsigned i = 0; i < set.count; ++i) {
      double a = 0.0, t = 0.0;
      for (unsigned v = 0; v < set.polygon[i].count; ++v) {
        a += set.polygon[i].vertex[v].inverse_velocity;
        t += set.polygon[i].vertex[v].time;
      }
      a /= (double) set.polygon[i].count;
      t /= (double) set.polygon[i].count;
      CHECK_EQUAL_INT(chopper_polygon_set_contains(&set, a, t), 1);
    }

    /* Over the whole sampled rectangle, containment and the independent point test must
     * agree everywhere except within rounding of the boundary. */
    unsigned disagreements = 0, inside = 0;
    for (unsigned i = 1; i <= 20000; ++i) {
      const double a = SOURCE_INVERSE_VELOCITY_MIN
                     + SOURCE_INVERSE_VELOCITY_RANGE * radical_inverse(i, 2);
      const double t = SOURCE_TIME_MIN + SOURCE_TIME_RANGE * radical_inverse(i, 3);
      const int by_polygon = chopper_polygon_set_contains(&set, a, t);
      const int by_disks = neutron_passes(bifrost, 6, a, t, 0.0);
      if (by_polygon) ++inside;
      if (by_polygon != by_disks) ++disagreements;
    }
    CHECK(inside > 0);
    CHECK_EQUAL_INT(disagreements, 0);
    chopper_polygon_set_free(&set);
  }

  TEST("the region writes itself as JSON");
  {
    chopper_polygon_set set = bifrost_source();
    const chopper_polygon sampled = chopper_polygon_rectangle(
        SOURCE_INVERSE_VELOCITY_MIN, SOURCE_INVERSE_VELOCITY_RANGE,
        SOURCE_TIME_MIN, SOURCE_TIME_RANGE);
    chopper_polygon_set_transmit_train(&set, 6, bifrost, NULL);

    CHECK_EQUAL_INT(chopper_write_polygons_to_file(".", "test_polygon_output", ".json",
                                                   "/", &set, &sampled), 1);
    FILE * written = fopen("./test_polygon_output.json", "r");
    CHECK(written != NULL);
    if (written) {
      /* Braces and brackets balanced, and the expected keys present: enough to catch a
       * malformed write without linking a JSON parser into the test. The Python side
       * parses it properly. */
      int braces = 0, brackets = 0, quotes = 0, saw_acceptance = 0, saw_vertices = 0;
      char buffer[4096];
      while (fgets(buffer, sizeof(buffer), written)) {
        for (const char * c = buffer; *c; ++c) {
          if (*c == '"') quotes ^= 1;
          if (quotes) continue;
          if (*c == '{') ++braces;
          if (*c == '}') --braces;
          if (*c == '[') ++brackets;
          if (*c == ']') --brackets;
        }
        if (strstr(buffer, "\"acceptance\"")) saw_acceptance = 1;
        if (strstr(buffer, "\"vertices\"")) saw_vertices = 1;
      }
      fclose(written);
      CHECK_EQUAL_INT(braces, 0);
      CHECK_EQUAL_INT(brackets, 0);
      CHECK_EQUAL_INT(saw_acceptance, 1);
      CHECK_EQUAL_INT(saw_vertices, 1);
      remove("./test_polygon_output.json");
    }

    /* No sampled region to measure against: written as null, not as a division by zero. */
    CHECK_EQUAL_INT(chopper_write_polygons_to_file(".", "test_polygon_null", ".json",
                                                   "/", &set, NULL), 1);
    remove("./test_polygon_null.json");
    chopper_polygon_set_free(&set);
  }

  TEST("an empty set gives a sampler with nothing in it");
  {
    chopper_polygon_set set = chopper_polygon_set_empty();
    chopper_polygon_sampler sampler = chopper_polygon_sampler_make(&set, 1.0);
    CHECK_EQUAL_INT(sampler.count, 0);
    CHECK_CLOSE(sampler.acceptance, 0.0, 0.0);
    CHECK(sampler.cumulative == NULL);
    double a = -1, t = -1;
    chopper_polygon_sampler_draw(&sampler, 0.5, 0.5, 0.5, &a, &t);   /* must not crash */
    chopper_polygon_sampler_free(&sampler);
    chopper_polygon_sampler_free(&sampler);   /* twice is fine */
    chopper_polygon_sampler_free(NULL);
  }

  return chopper_test_report();
}
