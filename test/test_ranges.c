/* The range algebra the chopper calculations are built out of.
 *
 * `chopper_inverse_velocity_windows` works by intersecting one set of allowed ranges per
 * chopper, so an error in the intersection is an error in every chopper answer.
 */
#include <stdlib.h>
#include "chopper-lib.h"
#include "test_util.h"

static void test_overlap_classification(void) {
  TEST("two ranges are classified by how they overlap");
  const range a = {0.0, 1.0};

  const range disjoint = {2.0, 3.0};
  CHECK_EQUAL_INT(classify_range_overlap(&a, &disjoint), 0);
  CHECK_EQUAL_INT(classify_range_overlap(&disjoint, &a), 0);

  const range same = {0.0, 1.0};
  CHECK_EQUAL_INT(classify_range_overlap(&a, &same), 1);

  const range inside = {0.25, 0.75};
  CHECK_EQUAL_INT(classify_range_overlap(&a, &inside), 3);  /* a contains b */
  CHECK_EQUAL_INT(classify_range_overlap(&inside, &a), 2);  /* b contains a */

  const range higher = {0.5, 1.5};
  CHECK_EQUAL_INT(classify_range_overlap(&a, &higher), -2); /* b reaches higher */
  CHECK_EQUAL_INT(classify_range_overlap(&higher, &a), -3); /* a reaches higher */
}

static void test_touching_ranges_count_as_overlapping(void) {
  TEST("ranges that meet at a point are not disjoint");
  const range a = {0.0, 1.0};
  const range b = {1.0, 2.0};
  CHECK(classify_range_overlap(&a, &b) != 0);
}

static void test_intersection_keeps_only_shared_spans(void) {
  TEST("an intersection keeps only what both sets cover");
  range a_ranges[2] = {{0.0, 1.0}, {2.0, 3.0}};
  range b_ranges[2] = {{0.5, 2.5}, {2.75, 4.0}};
  range_set a = {2, a_ranges};
  range_set b = {2, b_ranges};

  range_set got = range_intersection(a, b);

  /* [0,1]&[0.5,2.5] -> [0.5,1]; [2,3]&[0.5,2.5] -> [2,2.5]; [2,3]&[2.75,4] -> [2.75,3] */
  CHECK_EQUAL_INT(got.count, 3);
  if (got.count == 3) {
    CHECK_CLOSE(got.ranges[0].minimum, 0.5, 1e-12);
    CHECK_CLOSE(got.ranges[0].maximum, 1.0, 1e-12);
    CHECK_CLOSE(got.ranges[1].minimum, 2.0, 1e-12);
    CHECK_CLOSE(got.ranges[1].maximum, 2.5, 1e-12);
    CHECK_CLOSE(got.ranges[2].minimum, 2.75, 1e-12);
    CHECK_CLOSE(got.ranges[2].maximum, 3.0, 1e-12);
  }
  if (got.ranges && got.ranges != a_ranges && got.ranges != b_ranges) free(got.ranges);
}

static void test_disjoint_sets_intersect_to_nothing(void) {
  TEST("sets that share no span intersect to nothing");
  range a_ranges[1] = {{0.0, 1.0}};
  range b_ranges[1] = {{2.0, 3.0}};
  range_set a = {1, a_ranges};
  range_set b = {1, b_ranges};

  range_set got = range_intersection(a, b);
  CHECK_EQUAL_INT(got.count, 0);
  if (got.ranges && got.ranges != a_ranges && got.ranges != b_ranges) free(got.ranges);
}

static void test_sorting_merges_what_it_can(void) {
  TEST("sorting a set orders it and merges what touches");
  range ranges[3] = {{2.0, 3.0}, {0.0, 1.0}, {0.5, 1.5}};
  range_set unsorted = {3, ranges};

  range_set sorted = range_set_sort(unsorted);

  /* [0,1] and [0.5,1.5] overlap and become [0,1.5]; [2,3] stands alone */
  CHECK_EQUAL_INT(sorted.count, 2);
  if (sorted.count == 2) {
    CHECK_CLOSE(sorted.ranges[0].minimum, 0.0, 1e-12);
    CHECK_CLOSE(sorted.ranges[0].maximum, 1.5, 1e-12);
    CHECK_CLOSE(sorted.ranges[1].minimum, 2.0, 1e-12);
    CHECK_CLOSE(sorted.ranges[1].maximum, 3.0, 1e-12);
  }
  if (sorted.ranges && sorted.ranges != ranges) free(sorted.ranges);
}

int main(void) {
  test_overlap_classification();
  test_touching_ranges_count_as_overlapping();
  test_intersection_keeps_only_shared_spans();
  test_disjoint_sets_intersect_to_nothing();
  test_sorting_merges_what_it_can();
  return chopper_test_report();
}
