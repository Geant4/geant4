# - Provide a custom target for validating sources.cmake and on-disk files
# As sources.cmake lists source files of Geant4 explicitely, we can often
# get a mismatch between this list and what's actually on disk.
# This module provides a custom target which executes a CMake script to check 
# for these errors and report mismatches. It fails with FATAL_ERROR if any
# mismatch is found, but will not do so until it has reported all errors.
#

# Configure the script
configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/Templates/geant4_validate_sources.cmake.in
  ${PROJECT_BINARY_DIR}/geant4_validate_sources.cmake
  @ONLY
  )

# Create the target
add_custom_target(validate_sources
  COMMAND ${CMAKE_COMMAND} -P ${PROJECT_BINARY_DIR}/geant4_validate_sources.cmake
  COMMENT "Validating Geant4 module sources"
  )



