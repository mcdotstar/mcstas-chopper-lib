# C tests, run by CTest. The Python tests alongside them build and run whole McStas
# instruments and are driven by pytest instead; these cover the library on its own.
foreach (name IN ITEMS test_ranges test_chopper_timing)
  add_executable(${name} ${name}.c)
  target_link_libraries(${name} PRIVATE chopper_lib)
  target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  add_test(NAME ${name} COMMAND ${name})
endforeach ()

if (WIN32)
  # The tests link chopper_lib as a shared library, so it has to be alongside the .exe
  foreach (name IN ITEMS test_ranges test_chopper_timing)
    add_custom_command(TARGET ${name} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              $<TARGET_FILE:chopper_lib> $<TARGET_FILE_DIR:${name}>)
  endforeach ()
endif ()
