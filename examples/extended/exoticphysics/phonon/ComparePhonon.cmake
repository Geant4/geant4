# - Compare test to reference file
#
# Run as:
#
# $ cmake -Dtest_file=<output> -Dreference_file=<ref> -P CompareFiles.cmake
#
# - Both files are checked for existence.
# - Comparison is via cmake -E compare_files
# - If comparison fails, print test file to stdout
#

# CL checks
set(__CF_USAGE "usage: cmake -Dtest_file=<input> -Dreference_file=<ref> -P CompareFiles.cmake")

if(NOT test_file)
  message(FATAL_ERROR "No test_file argument\n${__CF_USAGE}")
endif()

if(NOT reference_file)
  message(FATAL_ERROR "No reference file argument\n${__CF_USAGE}")
endif()

# File existence checks
# Assume relative paths are relative to cwd
if(NOT (IS_ABSOLUTE "${test_file}"))
  set(test_file "${CMAKE_CURRENT_SOURCE_DIR}/${test_file}")
endif()

if(NOT (IS_ABSOLUTE "${reference_file}"))
  set(reference_file "${CMAKE_CURRENT_SOURCE_DIR}/${reference_file}")
endif()

if(NOT EXISTS "${test_file}")
  message(FATAL_ERROR "test_file '${test_file}' does not exist")
endif()

if(NOT (EXISTS "${reference_file}"))
  message(FATAL_ERROR "reference_file '${reference_file}' does not exist")
endif()

# File comparison
execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files ${test_file} ${reference_file}
  RESULT_VARIABLE __cf_result_var
  OUTPUT_QUIET
  ERROR_QUIET
  )

if(__cf_result_var)
  # Comparison failed, print contents of test_file and exit with an error
  file(READ ${test_file} __cf_test_contents)
  message("-- Test output '${test_file}' does not match reference '${reference_file}'")
  message("-- Contents of test output '${test_file}' (unformatted):")
  message("${__cf_test_contents}")
  message(FATAL_ERROR "TEST FAILED")
else()
  message(STATUS "Test output '${test_file}' matches reference '${reference_file}'")
endif()

