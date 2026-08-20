#ifndef CHOPPER_LIB_TEST_UTIL_H
#define CHOPPER_LIB_TEST_UTIL_H
/* A test harness small enough to not be a dependency.
 *
 * Every check reports where it failed and keeps going, so one run tells you everything
 * that is broken rather than only the first thing.
 */
#include <stdio.h>
#include <math.h>

static int chopper_test_failures = 0;
static const char * chopper_test_name = "";

#define TEST(name) \
  chopper_test_name = (name); \
  printf("-- %s\n", chopper_test_name);

#define CHECK(condition) \
  do { \
    if (!(condition)) { \
      ++chopper_test_failures; \
      printf("   FAIL %s:%d in %s: %s\n", __FILE__, __LINE__, chopper_test_name, #condition); \
    } \
  } while (0)

#define CHECK_CLOSE(actual, expected, tolerance) \
  do { \
    const double a_ = (actual), e_ = (expected), t_ = (tolerance); \
    if (!(fabs(a_ - e_) <= t_)) { \
      ++chopper_test_failures; \
      printf("   FAIL %s:%d in %s: %s is %.12g, expected %.12g (tolerance %.3g)\n", \
             __FILE__, __LINE__, chopper_test_name, #actual, a_, e_, t_); \
    } \
  } while (0)

#define CHECK_EQUAL_INT(actual, expected) \
  do { \
    const long a_ = (long)(actual), e_ = (long)(expected); \
    if (a_ != e_) { \
      ++chopper_test_failures; \
      printf("   FAIL %s:%d in %s: %s is %ld, expected %ld\n", \
             __FILE__, __LINE__, chopper_test_name, #actual, a_, e_); \
    } \
  } while (0)

static int chopper_test_report(void) {
  if (chopper_test_failures) {
    printf("\n%d check(s) failed\n", chopper_test_failures);
    return 1;
  }
  printf("\nall checks passed\n");
  return 0;
}

#endif //CHOPPER_LIB_TEST_UTIL_H
