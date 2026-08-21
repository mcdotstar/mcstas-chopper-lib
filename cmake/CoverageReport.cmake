# Turn the counters CTest just produced into a report. Run in script mode by the
# `coverage` target; not meant to be included.
#
# Expects: BINARY_DIR, SOURCE_DIR, OUTPUT_DIR, GCOV_TOOL, and optionally GCOVR/LCOV/GENHTML.
#
# gcovr and lcov each produce a browsable HTML report, but neither is guaranteed to be
# installed. gcov is, because it ships with the compiler, so it is the fallback and the
# report degrades to a per-file text summary rather than disappearing.

file(GLOB_RECURSE gcda_files "${BINARY_DIR}/*.gcda")
if (NOT gcda_files)
  message(FATAL_ERROR
    "No .gcda files under ${BINARY_DIR}.\n"
    "Coverage counters are written when an instrumented test runs, so either the build "
    "was not configured with -DCHOPPER_LIB_COVERAGE=ON or no test ran.")
endif ()

file(MAKE_DIRECTORY "${OUTPUT_DIR}")

if (GCOVR)
  message(STATUS "Coverage report via gcovr")
  execute_process(
    COMMAND ${GCOVR}
            --root "${SOURCE_DIR}"
            --gcov-executable "${GCOV_TOOL}"
            # the library is what is under test; the harness that exercises it is not
            --filter "${SOURCE_DIR}/chopper-lib\\.c"
            --print-summary
            --txt "${OUTPUT_DIR}/coverage.txt"
            --html-details "${OUTPUT_DIR}/index.html"
            --xml "${OUTPUT_DIR}/coverage.xml"
    WORKING_DIRECTORY "${BINARY_DIR}"
    RESULT_VARIABLE status)
  if (NOT status EQUAL 0)
    message(FATAL_ERROR "gcovr failed (${status})")
  endif ()
  message(STATUS "HTML report: ${OUTPUT_DIR}/index.html")

elseif (LCOV AND GENHTML)
  message(STATUS "Coverage report via lcov")
  execute_process(
    COMMAND ${LCOV} --capture --directory "${BINARY_DIR}" --base-directory "${SOURCE_DIR}"
                    --gcov-tool "${GCOV_TOOL}" --output-file "${OUTPUT_DIR}/coverage.info"
                    --rc branch_coverage=1 --quiet
    RESULT_VARIABLE status)
  if (NOT status EQUAL 0)
    message(FATAL_ERROR "lcov --capture failed (${status})")
  endif ()
  execute_process(
    COMMAND ${LCOV} --extract "${OUTPUT_DIR}/coverage.info" "${SOURCE_DIR}/chopper-lib.c"
                    --output-file "${OUTPUT_DIR}/coverage.info" --quiet
    RESULT_VARIABLE status)
  if (NOT status EQUAL 0)
    message(FATAL_ERROR "lcov --extract failed (${status})")
  endif ()
  execute_process(COMMAND ${LCOV} --list "${OUTPUT_DIR}/coverage.info")
  execute_process(
    COMMAND ${GENHTML} "${OUTPUT_DIR}/coverage.info" --output-directory "${OUTPUT_DIR}" --quiet)
  message(STATUS "HTML report: ${OUTPUT_DIR}/index.html")

else ()
  message(STATUS "Coverage report via gcov (install gcovr or lcov for an HTML report)")
  # gcov writes its .gcov files into the working directory, so give it one of its own
  file(REMOVE_RECURSE "${OUTPUT_DIR}/gcov")
  file(MAKE_DIRECTORY "${OUTPUT_DIR}/gcov")
  set(summary "")
  foreach (gcda IN LISTS gcda_files)
    get_filename_component(gcda_dir "${gcda}" DIRECTORY)
    execute_process(
      COMMAND ${GCOV_TOOL} --branch-probabilities --object-directory "${gcda_dir}" "${gcda}"
      WORKING_DIRECTORY "${OUTPUT_DIR}/gcov"
      OUTPUT_VARIABLE out ERROR_VARIABLE err RESULT_VARIABLE status)
    if (NOT status EQUAL 0)
      message(FATAL_ERROR "gcov failed on ${gcda} (${status}):\n${err}")
    endif ()
    # "File 'x.c'" followed by "Lines executed:12.34% of 567"
    string(REGEX MATCHALL "File '[^']*'\n[^\n]*Lines executed:[0-9.]+% of [0-9]+" hits "${out}")
    foreach (hit IN LISTS hits)
      string(REGEX MATCH "File '([^']*)'" _ "${hit}")
      set(file "${CMAKE_MATCH_1}")
      string(REGEX MATCH "Lines executed:([0-9.]+)% of ([0-9]+)" _ "${hit}")
      # only the library's own sources; not the test harness or system headers
      get_filename_component(base "${file}" NAME)
      if (base STREQUAL "chopper-lib.c" OR base STREQUAL "chopper-lib.h")
        list(APPEND summary "  ${base}: ${CMAKE_MATCH_1}% of ${CMAKE_MATCH_2} lines")
      endif ()
    endforeach ()
  endforeach ()
  if (NOT summary)
    message(FATAL_ERROR "gcov produced no data for chopper-lib.c")
  endif ()
  list(REMOVE_DUPLICATES summary)
  message(STATUS "Line coverage:")
  foreach (line IN LISTS summary)
    message(STATUS "${line}")
  endforeach ()
  message(STATUS "Annotated sources: ${OUTPUT_DIR}/gcov")
endif ()
