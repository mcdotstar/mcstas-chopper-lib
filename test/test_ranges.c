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

/* An empty set sorts to itself, and nothing is dereferenced on the way.
 *
 * `{0, NULL}` is ordinary here, not exceptional: it is what a chopper admitting no
 * window intersects with, and what a train that passes nothing carries from there on.
 * Sorting it used to reach `qsort` with a null pointer, which is undefined however
 * little there is to sort -- run this suite under UBSan to see the difference.
 */
static void test_sorting_an_empty_set_is_a_no_op(void) {
  TEST("an empty set sorts to itself");
  range_set empty = {0, NULL};

  range_set sorted = range_set_sort(empty);

  CHECK_EQUAL_INT(sorted.count, 0);
  CHECK(sorted.ranges == NULL);

  /* and it still intersects with a real set, to nothing */
  range ranges[1] = {{0.0, 1.0}};
  range_set real = {1, ranges};
  range_set got = range_intersection(empty, real);
  CHECK_EQUAL_INT(got.count, 0);
  if (got.ranges && got.ranges != ranges) free(got.ranges);
}

static void test_a_contained_range_does_not_shrink_the_one_it_is_in(void) {
  TEST("merging keeps the wider range, whichever order the two arrive in");
  /* This was [1, 3] both ways round: the merge assigned the later range's upper edge to
   * the range being built rather than taking the larger of the two. */
  range wide_first[2] = {{1.0, 10.0}, {2.0, 3.0}};
  range_set a = range_set_sort((range_set){2, wide_first});
  CHECK_EQUAL_INT(a.count, 1);
  CHECK_CLOSE(a.ranges[0].minimum, 1.0, 0.0);
  CHECK_CLOSE(a.ranges[0].maximum, 10.0, 0.0);

  range narrow_first[2] = {{2.0, 3.0}, {1.0, 10.0}};
  range_set b = range_set_sort((range_set){2, narrow_first});
  CHECK_EQUAL_INT(b.count, 1);
  CHECK_CLOSE(b.ranges[0].minimum, 1.0, 0.0);
  CHECK_CLOSE(b.ranges[0].maximum, 10.0, 0.0);
}

/* Every ordering of the same set must merge to the same answer.
 *
 * qsort is not stable and orders ties as each platform's implementation pleases, so a set
 * with tied lower edges used to merge differently on Windows and on Linux -- which is how
 * a chopper train's admitted band came out wrong on one of them. Permuting the input here
 * stands in for that, without needing the other platform to run it on.
 */
static void permute_and_check(range * work, range * source, unsigned n, unsigned depth,
                              unsigned * used, range * expected, unsigned expected_count) {
  if (depth == n) {
    range copy[8];
    for (unsigned i = 0; i < n; ++i) copy[i] = work[i];
    range_set got = range_set_sort((range_set){n, copy});
    CHECK_EQUAL_INT(got.count, expected_count);
    if (got.count == expected_count) {
      for (unsigned i = 0; i < expected_count; ++i) {
        CHECK_CLOSE(got.ranges[i].minimum, expected[i].minimum, 0.0);
        CHECK_CLOSE(got.ranges[i].maximum, expected[i].maximum, 0.0);
      }
    }
    return;
  }
  for (unsigned i = 0; i < n; ++i) {
    if (used[i]) continue;
    used[i] = 1;
    work[depth] = source[i];
    permute_and_check(work, source, n, depth + 1, used, expected, expected_count);
    used[i] = 0;
  }
}

static void test_the_merged_answer_does_not_depend_on_the_input_order(void) {
  TEST("every ordering of a set with tied and nested ranges merges the same way");
  /* Tied lower edges and a contained range, which is what clamping to the search bounds
   * produces in chopper_inverse_velocity_windows. */
  range source[5] = {{0.0, 10.0}, {0.0, 1.0}, {5.0, 6.0}, {20.0, 25.0}, {24.0, 24.5}};
  range expected[2] = {{0.0, 10.0}, {20.0, 25.0}};
  range work[5];
  unsigned used[5] = {0, 0, 0, 0, 0};
  permute_and_check(work, source, 5, 0, used, expected, 2);
}

static void test_touching_and_disjoint_ranges_still_behave(void) {
  TEST("touching ranges merge and separated ones do not");
  range touching[2] = {{0.0, 1.0}, {1.0, 2.0}};
  range_set merged = range_set_sort((range_set){2, touching});
  CHECK_EQUAL_INT(merged.count, 1);
  CHECK_CLOSE(merged.ranges[0].maximum, 2.0, 0.0);

  range apart[2] = {{0.0, 1.0}, {1.5, 2.0}};
  range_set kept = range_set_sort((range_set){2, apart});
  CHECK_EQUAL_INT(kept.count, 2);
}

static void test_an_inverted_range_is_put_the_right_way_round(void) {
  TEST("a range given upper edge first is normalised before merging");
  /* The set sort normalises its own ranges. There used to be a `range_sort` for this
   * which took its argument by value, and so never did it; 4.2.1 removed it. */
  range inverted[2] = {{10.0, 1.0}, {2.0, 3.0}};
  range_set got = range_set_sort((range_set){2, inverted});
  CHECK_EQUAL_INT(got.count, 1);
  CHECK_CLOSE(got.ranges[0].minimum, 1.0, 0.0);
  CHECK_CLOSE(got.ranges[0].maximum, 10.0, 0.0);
}

int main(void) {
  test_overlap_classification();
  test_touching_ranges_count_as_overlapping();
  test_intersection_keeps_only_shared_spans();
  test_disjoint_sets_intersect_to_nothing();
  test_sorting_merges_what_it_can();
  test_sorting_an_empty_set_is_a_no_op();
  test_a_contained_range_does_not_shrink_the_one_it_is_in();
  test_the_merged_answer_does_not_depend_on_the_input_order();
  test_touching_and_disjoint_ranges_still_behave();
  test_an_inverted_range_is_put_the_right_way_round();
  return chopper_test_report();
}
