# - Script for configuring and installing a Modulefile for Geant4
#
# Environment Modules is a standard tool for configuring the environment
# for a package in a shel/intepreter agnostic way. See:
#
#  http://modules.sourceforge.net/
#
# As with other Geant4 tool support, a template file is provided to
# generate the modulefile using the known build settings. Because
# modulefiles are generally only used on installed packages, only
# an Install Tree file is generated.
#
# The resultant modulefile is installed to the share directory, though
# this is not intended to be its final location. System Admins may
# wish to move the file to their local modulefile path
#

macro(_geant4_prepare_modulefile_inputs)
  # - Compatibility
  if(APPLE)
    set(DYNAMIC_LOADER_PATHNAME "DYLD_LIBRARY_PATH")
  else()
    set(DYNAMIC_LOADER_PATHNAME "LD_LIBRARY_PATH")
  endif()

  # - Datasets
  geant4_get_datasetnames(GEANT4_EXPORTED_DATASETS)
  set(G4DATASET_TCLLIST)
  list(REMOVE_ITEM GEANT4_EXPORTED_DATASETS "G4ENSDFSTATE")
  foreach(_ds ${GEANT4_EXPORTED_DATASETS})
    geant4_get_dataset_property(${_ds} ENVVAR ${_ds}_ENVVAR)
    geant4_get_dataset_property(${_ds} INSTALL_DIR ${_ds}_PATH)
    set(G4DATASET_TCLLIST "${G4DATASET_TCLLIST} ${${_ds}_ENVVAR} ${${_ds}_PATH}")
  endforeach()
endmacro()


function(geant4_configure_modulefile)
  # - prepare configuration environment
  _geant4_prepare_modulefile_inputs()

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
endfunction()

