# - Setup for installing architecture independent read only files.
#
# There are only two main items to install here:
#
#  Geant4 Examples - basically everything under the 'examples' directory.
#
#  Geant4 Data Libraries - not supplied with code, these must be downloaded
#                          or installed from a local URL.
#


#----------------------------------------------------------------------------
# Install examples if requested
#
option(GEANT4_INSTALL_EXAMPLES "Install source code for Geant4 examples" OFF)
mark_as_advanced(GEANT4_INSTALL_EXAMPLES)

if(GEANT4_INSTALL_EXAMPLES)
  install(DIRECTORY examples
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Geant4-${Geant4_VERSION}
    COMPONENT Examples
    PATTERN "CVS" EXCLUDE
    PATTERN ".svn" EXCLUDE)
endif()

GEANT4_ADD_FEATURE(GEANT4_INSTALL_EXAMPLES "Will install source code for Geant4 examples")

#----------------------------------------------------------------------------
# Install data libraries if requested
# Because we use ExternalProject to drive the download and unpacking, 
# we restrict the availability of the option to CMake >= 2.8
# We check by requiring CMAKE_VERSION to be greater than 2.7. Since a 2.7
# version of CMake has never been released, this should be o.k.!
# If we don't have the required version, warn and unset the variable.
# This handles the case where users may pass the variable on the command line
if(NOT ${CMAKE_VERSION} VERSION_GREATER 2.7)
  if(GEANT4_INSTALL_DATA)
    message(STATUS "WARNING: GEANT4_INSTALL_DATA requires CMake >= 2.8 -- deactivating")
  endif()
  set(GEANT4_INSTALL_DATA OFF CACHE BOOL "Install data NOT supported on CMake <2.8" FORCE)
else()
  option(GEANT4_INSTALL_DATA "Download and install Geant4 Data Libraries" OFF)

  # For very large datasets, and on CMake 2.8.1 and above only, provide
  # a cached variable to allow the TIMEOUT to be adjusted if needed for slow
  # connections. Default to 1500 seconds (25 minutes)
  if(${CMAKE_VERSION} VERSION_GREATER 2.8.0)
    set(GEANT4_INSTALL_DATA_TIMEOUT 1500 CACHE STRING "Time in seconds allowed for Data Library download")
    set(_g4id_timeout_cmd TIMEOUT ${GEANT4_INSTALL_DATA_TIMEOUT})
  endif()


  if(GEANT4_INSTALL_DATA)
    include(ExternalProject)
    set(_urlprefix "http://geant4.cern.ch/support/source")
    set(GEANT4_DATASETS
      G4NDL/4.0/G4NDL/tar.gz/G4NEUTRONHPDATA
      G4EMLOW/6.26/G4EMLOW/tar.gz/G4LEDATA
      PhotonEvaporation/2.2/PhotonEvaporation/tar.gz/G4LEVELGAMMADATA
      RadioactiveDecay/3.4/G4RadioactiveDecay/tar.gz/G4RADIOACTIVEDATA
      G4ABLA/3.0/G4ABLA/tar.gz/G4ABLADATA
      G4NEUTRONXS/1.1/G4NEUTRONXS/tar.gz/G4NEUTRONXSDATA
      G4PII/1.3/G4PII/tar.gz/G4PIIDATA
      RealSurface/1.0/RealSurface/tar.gz/G4REALSURFACEDATA
    )

    foreach(_ds ${GEANT4_DATASETS})
      string(REPLACE "/" ";" _tuple ${_ds})
      list(GET _tuple 0 _name)
      list(GET _tuple 1 _vers)
      list(GET _tuple 2 _fnam)
      list(GET _tuple 3 _suffix)

      ExternalProject_Add(${_name} 
        URL ${_urlprefix}/${_fnam}.${_vers}.${_suffix}
        ${_g4id_timeout_cmd}
        PREFIX Externals/${_fnam}-${_vers}
        SOURCE_DIR data/${_name}${_vers}
        CONFIGURE_COMMAND "" 
        BUILD_COMMAND "" 
        INSTALL_COMMAND ""
      )

      install(DIRECTORY ${CMAKE_BINARY_DIR}/data/${_name}${_vers}
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Geant4-${Geant4_VERSION}/data
        COMPONENT Data
      )

    endforeach()
  endif()

  GEANT4_ADD_FEATURE(GEANT4_INSTALL_DATA "Will download and install data libraries")
endif()




