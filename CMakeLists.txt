cmake_minimum_required(VERSION 3.25)
project(chopper_lib C)

set(CMAKE_C_STANDARD 11)

option(CHOPPER_LIB_COVERAGE "Instrument the library and its tests to measure line coverage" OFF)

# The unit conversions chopper-lib.c uses, which it does not define itself: McStas
# defines all three in its runtime, and the source is copied verbatim into every
# instrument that %includes it, so a definition here would be a second definition of
# PI in the generated C. Compiling the source anywhere else has to supply them, and
# the file refuses to compile without them rather than fall back to its own numbers.
# The values are the McStas runtime's. Cached, so that a consumer pulling this in with
# FetchContent can read them: one that compiles chopper-lib.c into a target of its own,
# instead of linking chopper_lib, needs
#   target_compile_definitions(its_target PRIVATE ${CHOPPER_LIB_DEFINITIONS})
set(CHOPPER_LIB_DEFINITIONS
    V2K=1.58825361e-3          # v [m/s] to k [1/AA]
    K2V=629.622368             # k [1/AA] to v [m/s]
    PI=3.14159265358979323846
    CACHE STRING "Constants chopper-lib.c is compiled with, matching the McStas runtime")

add_library(chopper_lib SHARED chopper-lib.c)
target_include_directories(chopper_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# PRIVATE: only the implementation reads them, and PI is too common a name to push onto
# everything that links the library.
target_compile_definitions(chopper_lib PRIVATE ${CHOPPER_LIB_DEFINITIONS})
# MSVC exports nothing from a DLL without __declspec(dllexport); the tests link against
# this library, so let CMake generate the export table rather than annotate every symbol.
set_target_properties(chopper_lib PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
if (NOT MSVC)
  # chopper-lib.c uses floor, ceil, fmod and fabs
  target_link_libraries(chopper_lib PUBLIC m)
endif ()

if (CHOPPER_LIB_COVERAGE)
  if (NOT CMAKE_C_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    message(FATAL_ERROR
      "CHOPPER_LIB_COVERAGE needs GCC or Clang; ${CMAKE_C_COMPILER_ID} has no --coverage.")
  endif ()
  # PUBLIC so the tests, which link this library, are instrumented too: the counters
  # measure which library lines the tests reached, and both halves have to agree.
  # Optimisation folds and reorders lines, which makes a line-level report misleading.
  target_compile_options(chopper_lib PUBLIC --coverage -O0 -g)
  target_link_options(chopper_lib PUBLIC --coverage)
endif ()

# Only when this is the project being built. Consumers pull the library in with
# FetchContent -- chopper-calculations does -- and should not inherit its test suite.
if (PROJECT_IS_TOP_LEVEL)
  include(CTest)
  if (BUILD_TESTING)
    add_subdirectory(test)

    if (CHOPPER_LIB_COVERAGE)
      # gcov ships with the compiler; the rest are nicer if they happen to be installed.
      if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
        find_program(GCOV_TOOL NAMES gcov-${CMAKE_C_COMPILER_VERSION} gcov REQUIRED)
      else ()
        find_program(LLVM_COV_TOOL NAMES llvm-cov REQUIRED)
        # gcovr and lcov both want a program that behaves like gcov; llvm-cov does only
        # under its `gcov` subcommand, so hand them a one-line shim that supplies it.
        set(GCOV_TOOL ${CMAKE_CURRENT_BINARY_DIR}/llvm-gcov.sh)
        file(WRITE ${GCOV_TOOL} "#!/bin/sh\nexec \"${LLVM_COV_TOOL}\" gcov \"$@\"\n")
        file(CHMOD ${GCOV_TOOL} PERMISSIONS
             OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE
             WORLD_READ WORLD_EXECUTE)
      endif ()
      find_program(GCOVR NAMES gcovr)
      find_program(LCOV NAMES lcov)
      find_program(GENHTML NAMES genhtml)

      set(coverage_dir ${CMAKE_CURRENT_BINARY_DIR}/coverage)
      add_custom_target(coverage
        # Counters accumulate across runs, so clear them or a stale run inflates the report
        COMMAND ${CMAKE_COMMAND} -E rm -rf ${coverage_dir}
        COMMAND ${CMAKE_COMMAND}
                -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/CoverageReset.cmake
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
        COMMAND ${CMAKE_COMMAND}
                -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
                -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
                -DOUTPUT_DIR=${coverage_dir}
                -DGCOV_TOOL=${GCOV_TOOL}
                -DGCOVR=${GCOVR}
                -DLCOV=${LCOV}
                -DGENHTML=${GENHTML}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/CoverageReport.cmake
        DEPENDS ${CHOPPER_LIB_TEST_TARGETS}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Running the tests and reporting line coverage"
        VERBATIM USES_TERMINAL)
    endif ()
  endif ()
endif ()
