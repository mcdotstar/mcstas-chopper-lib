# Delete the counters from any previous run. Run in script mode by the `coverage` target.
#
# gcov counters accumulate into .gcda across runs, so without this a report describes
# every test run since the build directory was created rather than this one.
file(GLOB_RECURSE gcda_files "${BINARY_DIR}/*.gcda")
if (gcda_files)
  file(REMOVE ${gcda_files})
endif ()
