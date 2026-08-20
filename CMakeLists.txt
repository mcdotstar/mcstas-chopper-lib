cmake_minimum_required(VERSION 3.25)
project(chopper_lib C)

set(CMAKE_C_STANDARD 11)

add_library(chopper_lib SHARED chopper-lib.c)
target_include_directories(chopper_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# MSVC exports nothing from a DLL without __declspec(dllexport); the tests link against
# this library, so let CMake generate the export table rather than annotate every symbol.
set_target_properties(chopper_lib PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
if (NOT MSVC)
  # chopper-lib.c uses floor, ceil, fmod and fabs
  target_link_libraries(chopper_lib PUBLIC m)
endif ()

# Only when this is the project being built. Consumers pull the library in with
# FetchContent -- chopper-calculations does -- and should not inherit its test suite.
if (PROJECT_IS_TOP_LEVEL)
  include(CTest)
  if (BUILD_TESTING)
    add_subdirectory(test)
  endif ()
endif ()
