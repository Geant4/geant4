# - Script for configuring and installing a Modulefile for Geant4
#
# Environment Modules is a standard tool for configuring the environment
# for a package in a shell/intepreter agnostic way. See:
#
#  http://modules.sourceforge.net/
#
# As with other Geant4 tool support, a template file is provided to
# generate the modulefile using the known build settings. Though
# modulefiles are generally only used for installed packages, modulefiles
# are generated for both the Build and Install Trees.
#
# The resultant modulefile for the Build Tree is only provided on an
# 'as is' basis. It is intended for Geant4 developers only, and is
# otherwise unsupported.
#
# The resultant modulefile for the Install Tree is installed to the
# share directory, though this is not intended to be its final location.
# System Admins may wish to move the file to their local modulefile path.
# Absolute paths to the install of Geant4 are used to help move the
# modulefile around, so if both modulefile and install of Geant4 are
# moved, the paths in the modulefile should be patched/sedded.
#

macro(_geant4_prepare_modulefile_inputs)
  cmake_parse_arguments(GPI
    ""
    "MODE"
    ""
    ${ARGN}
    )

  # - Paths
  if("${GPI_MODE}" STREQUAL "INSTALL")
    set(GEANT4_MODULEFILE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")
    set(GEANT4_MODULEFILE_INSTALL_BINDIR "${CMAKE_INSTALL_FULL_BINDIR}")
    set(GEANT4_MODULEFILE_INSTALL_LIBDIR "${CMAKE_INSTALL_FULL_LIBDIR}")
    geant4_export_datasets(INSTALL GEANT4_EXPORTED_DATASETS)
  else()
    set(GEANT4_MODULEFILE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}")
    set(GEANT4_MODULEFILE_INSTALL_BINDIR "${PROJECT_BINARY_DIR}")
    set(GEANT4_MODULEFILE_INSTALL_LIBDIR "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
    geant4_export_datasets(BUILD GEANT4_EXPORTED_DATASETS)
  endif()

  # - Compatibility
  if(APPLE)
    set(DYNAMIC_LOADER_PATHNAME "DYLD_LIBRARY_PATH")
  else()
    set(DYNAMIC_LOADER_PATHNAME "LD_LIBRARY_PATH")
  endif()

  # - Datasets
  set(G4DATASET_TCLLIST)
  foreach(_ds ${GEANT4_EXPORTED_DATASETS})
    # listify tuple
    string(REPLACE "|" ";" _ds "${_ds}")
    # Extract envar and path entries
    list(GET _ds 1 _ds_ENVVAR)
    list(GET _ds 2 _ds_PATH)
    set(G4DATASET_TCLLIST "${G4DATASET_TCLLIST} ${_ds_ENVVAR} ${_ds_PATH}")
  endforeach()
endmacro()


function(geant4_configure_modulefile)
  # Install Tree
  # - prepare configuration environment
  _geant4_prepare_modulefile_inputs(MODE INSTALL)

  # - Configure the file
  configure_file(
    ${PROJECT_SOURCE_DIR}/cmake/Templates/geant4-modulefile.in
    ${PROJECT_BINARY_DIR}/InstallTreeFiles/geant4-${Geant4_VERSION}
    @ONLY
    )

  # - Install it
  install(FILES  ${PROJECT_BINARY_DIR}/InstallTreeFiles/geant4-${Geant4_VERSION}
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Geant4-${Geant4_VERSION}
    COMPONENT Development
    )

  # Build Tree
  # - prepare configuration environment
  _geant4_prepare_modulefile_inputs()

  # - Configure the file
  configure_file(
    ${PROJECT_SOURCE_DIR}/cmake/Templates/geant4-modulefile.in
    ${PROJECT_BINARY_DIR}/geant4-${Geant4_VERSION}
    @ONLY
    )
endfunction()

