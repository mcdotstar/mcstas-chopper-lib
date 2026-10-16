/* How likely a neutron is to get through.
 *
 * Held to the definition rather than to the code: a neutron at position x across the beam
 * crosses each disk of a group at that same x, sees an edge `sweep * h * x` early or late,
 * and is passed by the group if every disk is open for it then. Averaging that over x on a
 * fine grid, and multiplying across groups, is the transmission -- worked out here from
 * the disk definitions directly, sharing nothing with chopper-lib's polygons.
 */
#include <stdlib.h>
#include <math.h>
#include "chopper-lib.h"
#include "test_util.h"

#define LATERAL 8000    /* positions across the beam in the reference: 2.5e-4 apart */

static double edge_time(const chopper_parameters d, const double a) {
  return d.delay + (d.beam - a) / 360.0 / d.speed;   /* as chopper-lib's header says */
}

static int open_at(const chopper_parameters d, const double arrival) {
  const double tau = 1.0 / fabs(d.speed);
  for (unsigned w = 0; w < d.edge_count / 2; ++w) {
    double o = edge_time(d, d.edges[2 * w]), c = edge_time(d, d.edges[2 * w + 1]);
    if (c < o) { const double s = o; o = c; c = s; }
    double phase = arrival - o;
    phase -= floor(phase / tau) * tau;
    if (phase <= c - o) return 1;
  }
  return 0;
}

static double reference_group(const chopper_parameters * d, const unsigned m,
                              const double u, const double t) {
  unsigned passed = 0;
  for (unsigned i = 0; i < LATERAL; ++i) {
    const double x = -1.0 + (i + 0.5) * 2.0 / LATERAL;
    int all = 1;
    for (unsigned k = 0; all && k < m; ++k) {
      const double h = d[k].aperture / 2.0 / 360.0 / fabs(d[k].speed);
      const double sweep = d[k].speed > 0 ? 1.0 : -1.0;
      all = open_at(d[k], t + d[k].path * u - sweep * h * x);
    }
    passed += all;
  }
  return (double) passed / LATERAL;
}

static double reference_train(const chopper_parameters * d, const unsigned n,
                              const unsigned * groups, const double u, const double t) {
  double product = 1.0;
  unsigned i = 0;
  while (i < n) {
    unsigned j = i + 1;
    while (groups && j < n && groups[j] == groups[i]) ++j;
    product *= reference_group(d + i, j - i, u, t);
    i = j;
  }
  return product;
}

static unsigned long long seed = 12345;
static double uniform(void) {
  seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
  return (double) (seed >> 11) / 9007199254740992.0;
}

static chopper_weighted_set rectangle(const double u0, const double du,
                                      const double t0, const double dt) {
  chopper_weighted_set set = chopper_weighted_set_empty();
  chopper_weighted_polygon piece;
  piece.polygon = chopper_polygon_rectangle(u0, du, t0, dt);
  piece.factor_count = 0;
  chopper_weighted_set_add(&set, &piece);
  return set;
}

/* The worst disagreement with the reference over points drawn across where the set is. */
static double worst_disagreement(const chopper_weighted_set * set,
                                 const chopper_parameters * d, const unsigned n,
                                 const unsigned * groups, const unsigned points) {
  double u_lo = 1e300, u_hi = -1e300, t_lo = 1e300, t_hi = -1e300;
  for (unsigned i = 0; i < set->count; ++i) {
    double lo, hi;
    chopper_polygon_extent(&set->polygon[i].polygon, 1.0, 0.0, &lo, &hi);
    if (lo < u_lo) u_lo = lo;
    if (hi > u_hi) u_hi = hi;
    chopper_polygon_extent(&set->polygon[i].polygon, 0.0, 1.0, &lo, &hi);
    if (lo < t_lo) t_lo = lo;
    if (hi > t_hi) t_hi = hi;
  }
  double worst = 0.0;
  for (unsigned p = 0; p < points; ++p) {
    const double u = u_lo + uniform() * (u_hi - u_lo);
    const double t = t_lo + uniform() * (t_hi - t_lo);
    const double diff = fabs(chopper_weighted_set_transmission(set, u, t)
                             - reference_train(d, n, groups, u, t));
    if (diff > worst) worst = diff;
  }
  return worst;
}

/* ---- the BIFROST train, as chopcal sets it for a 3 AA band, beam apertures and all ---- */
static double ps_edges[2] = {-85.0, 85.0};
static double fo1_edges[2] = {-19.13, 19.13};
static double fo2_edges[2] = {-26.005, 26.005};
static double bw_edges[2] = {-80.5, 80.5};
static chopper_parameters bifrost[6] = {
  {196.0, 0.0038781004374078398, 0.0, 2, ps_edges,  6.3420000000000005,  5.650725436158406},
  {196.0, 0.0060873974895620350, 0.0, 2, ps_edges,  6.3620000000000000,  5.650725436158406},
  { 14.0, 0.0061399701448323090, 0.0, 2, fo1_edges,  8.5300000000000000,  5.650725436158406},
  { 14.0, 0.0095476374535256970, 0.0, 2, fo2_edges, 14.9730000000000000,  5.650725436158406},
  { 14.0, 0.0428822715471184300, 0.0, 2, bw_edges, 78.0000000000000000, 13.379640178419610},
  {-14.0, 0.0428822715471184300, 0.0, 2, bw_edges, 78.0200000000000000, 13.379640178419610},
};
static const unsigned bifrost_groups[6] = {0, 0, 1, 2, 3, 3};

int main(void) {
  TEST("a product of affine factors integrates exactly");
  {
    chopper_weighted_polygon square;
    square.polygon = chopper_polygon_rectangle(0.0, 1.0, 0.0, 1.0);
    square.factor_count = 0;
    CHECK_CLOSE(chopper_weighted_polygon_integral(&square), 1.0, 1e-15);
    const chopper_affine x = {1.0, 0.0, 0.0}, y = {0.0, 1.0, 0.0}, one_minus_x = {-1.0, 0.0, 1.0};
    square.factor_count = 2; square.factor[0] = x; square.factor[1] = y;
    CHECK_CLOSE(chopper_weighted_polygon_integral(&square), 0.25, 1e-15);
    square.factor[1] = x;
    CHECK_CLOSE(chopper_weighted_polygon_integral(&square), 1.0 / 3.0, 1e-15);
    square.factor_count = 3; square.factor[2] = one_minus_x;   /* x^2 (1 - x) */
    CHECK_CLOSE(chopper_weighted_polygon_integral(&square), 1.0 / 12.0, 1e-15);
  }

  double lone_edges[2] = {-10.0, 10.0};
  const chopper_parameters lone = {14.0, 0.02, 0.0, 2, lone_edges, 0.0, 5.0};
  const double lone_width = 20.0 / 360.0 / 14.0;

  TEST("a lone disk's ramps integrate to its point-beam opening");
  {
    chopper_weighted_set set = rectangle(1e-4, 1e-4, 0.0, 0.05);
    CHECK_EQUAL_INT(chopper_weighted_set_transmit_group(&set, 1, &lone), 1);
    CHECK_CLOSE(chopper_weighted_set_integral(&set) / (lone_width * 1e-4), 1.0, 1e-12);
    chopper_weighted_set_free(&set);
  }

  TEST("an opening narrower than the beam never opens fully, and integrates the same");
  {
    double narrow_edges[2] = {-1.0, 1.0};
    const chopper_parameters narrow = {14.0, 0.02, 0.0, 2, narrow_edges, 0.0, 5.0};
    chopper_weighted_set set = rectangle(1e-4, 1e-4, 0.0, 0.05);
    CHECK_EQUAL_INT(chopper_weighted_set_transmit_group(&set, 1, &narrow), 1);
    CHECK_CLOSE(chopper_weighted_set_integral(&set) / (2.0 / 360.0 / 14.0 * 1e-4), 1.0,
                1e-12);
    CHECK_CLOSE(chopper_weighted_set_transmission(&set, 1.5e-4, 0.02), 0.4, 1e-12);
    chopper_weighted_set_free(&set);
  }

  TEST("a lone disk, either way round, is its definition");
  {
    chopper_parameters backwards = lone;
    backwards.speed = -14.0;
    const chopper_parameters both[2] = {lone, backwards};
    for (unsigned k = 0; k < 2; ++k) {
      chopper_weighted_set set = rectangle(1e-4, 1e-3, 0.0, 0.05);
      CHECK_EQUAL_INT(chopper_weighted_set_transmit_group(&set, 1, &both[k]), 1);
      CHECK(worst_disagreement(&set, &both[k], 1, NULL, 500) < 1e-3);
      chopper_weighted_set_free(&set);
    }
  }

  /* A pair sharing a housing, both with the opening centred on the beam at the same time.
   * Where their edges cross the beam centre, s = o, a co-rotating pair is half open -- the
   * smaller of two equal ramps -- and a counter-rotating pair is shut: its edges come in
   * from opposite sides and only meet there. */
  chopper_parameters pair_same[2] = {lone, lone};
  chopper_parameters pair_opposite[2] = {lone, lone};
  pair_opposite[1].speed = -14.0;
  const double o = edge_time(lone, 10.0), h = 5.0 / 2.0 / 360.0 / 14.0;

  TEST("a co-rotating pair passes the smaller of its two transmissions");
  {
    chopper_weighted_set set = rectangle(1e-4, 1e-4, 0.0, 0.05);
    CHECK_EQUAL_INT(chopper_weighted_set_transmit_group(&set, 2, pair_same), 1);
    CHECK_CLOSE(chopper_weighted_set_transmission(&set, 1.5e-4, o), 0.5, 1e-9);
    CHECK_CLOSE(chopper_weighted_set_transmission(&set, 1.5e-4, o + h / 2.0), 0.75, 1e-9);
    CHECK(worst_disagreement(&set, pair_same, 2, (const unsigned[]){0, 0}, 500) < 1e-3);
    chopper_weighted_set_free(&set);
  }

  TEST("a counter-rotating pair opens from the half-way time at twice the speed");
  {
    chopper_weighted_set set = rectangle(1e-4, 1e-4, 0.0, 0.05);
    CHECK_EQUAL_INT(chopper_weighted_set_transmit_group(&set, 2, pair_opposite), 1);
    CHECK_CLOSE(chopper_weighted_set_transmission(&set, 1.5e-4, o - h / 2.0), 0.0, 1e-9);
    CHECK_CLOSE(chopper_weighted_set_transmission(&set, 1.5e-4, o + h / 2.0), 0.5, 1e-9);
    CHECK_CLOSE(chopper_weighted_set_transmission(&set, 1.5e-4, o + h), 1.0, 1e-9);
    CHECK(worst_disagreement(&set, pair_opposite, 2, (const unsigned[]){0, 0}, 500) < 1e-3);
    chopper_weighted_set_free(&set);
  }

  TEST("the same pair as two independent groups multiplies instead");
  {
    chopper_weighted_set set = rectangle(1e-4, 1e-4, 0.0, 0.05);
    CHECK_EQUAL_INT(chopper_weighted_set_transmit_train(&set, 2, pair_opposite,
                                                        (const unsigned[]){0, 1}), 1);
    CHECK_CLOSE(chopper_weighted_set_transmission(&set, 1.5e-4, o), 0.25, 1e-9);
    CHECK(worst_disagreement(&set, pair_opposite, 2, (const unsigned[]){0, 1}, 500) < 1e-3);
    chopper_weighted_set_free(&set);
  }

  TEST("BIFROST, pairs sharing housings and the rest independent, is its definition");
  {
    chopper_weighted_set set = rectangle(1.0e-4, 1.13e-2, 0.0, 3.0e-3);
    CHECK_EQUAL_INT(chopper_weighted_set_transmit_train(&set, 6, bifrost, bifrost_groups), 1);
    CHECK(set.count > 0);
    CHECK(worst_disagreement(&set, bifrost, 6, bifrost_groups, 800) < 1.5e-3);

    /* never more than the support the same apertures give */
    chopper_polygon_set support = chopper_polygon_set_empty();
    const chopper_polygon r = chopper_polygon_rectangle(1.0e-4, 1.13e-2, 0.0, 3.0e-3);
    chopper_polygon_set_add(&support, &r);
    chopper_polygon_set_transmit_train(&support, 6, bifrost, NULL);
    const double integral = chopper_weighted_set_integral(&set);
    CHECK(integral > 0.0);
    CHECK(integral < chopper_polygon_set_area(&support));
    printf("   %u weighted pieces; integral %.6g of support %.6g s^2/m\n", set.count,
           integral, chopper_polygon_set_area(&support));
    chopper_polygon_set_free(&support);
    chopper_weighted_set_free(&set);
  }

  TEST("a group's disks must be adjacent");
  {
    chopper_weighted_set set = rectangle(1.0e-4, 1.13e-2, 0.0, 3.0e-3);
    CHECK_EQUAL_INT(chopper_weighted_set_transmit_train(&set, 3, bifrost,
                                                        (const unsigned[]){0, 1, 0}), 0);
    CHECK_EQUAL_INT(set.count, 0);
  }

  return chopper_test_report();
}
