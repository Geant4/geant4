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
option(GEANT4_INSTALL_EXAMPLES "Install all Geant4 examples" OFF)

if(GEANT4_INSTALL_EXAMPLES)
    install(DIRECTORY examples
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Geant4-${Geant4_VERSION}
        COMPONENT Examples
        PATTERN "CVS" EXCLUDE
        PATTERN ".svn" EXCLUDE)
endif()


#----------------------------------------------------------------------------
# Install data libraries if requested
# Because we use ExternalProject to drive the download and unpacking, 
# we restrict the availability of the option to CMake >= 2.8
# We check by requiring CMAKE_VERSION to be greater than 2.7. Since a 2.7
# version of CMake has never been released, this should be o.k.!
if(${CMAKE_VERSION} VERSION_GREATER 2.7)
  option(GEANT4_INSTALL_DATA "Download and install Geant4 Data Libraries" OFF)

  if(GEANT4_INSTALL_DATA)
    include(ExternalProject)
    set(_urlprefix "http://geant4.cern.ch/support/source")
    set(_datasets 
      G4NDL/3.16/G4NDL/tar.gz 
      G4EMLOW/6.23/G4EMLOW/tar.gz 
      PhotonEvaporation/2.2/G4PhotonEvaporation/tar.gz
      RadioactiveDecay/3.3/G4RadioactiveDecay/tar.gz
      G4ABLA/3.0/G4ABLA/tar.gz
      G4NEUTRONXS/1.1/G4NEUTRONXS/tar.gz
      G4PII/1.3/G4PII/tar.gz
      RealSurface/1.0/RealSurface/tar.gz
    )

    foreach(_ds ${_datasets})
      string(REPLACE "/" ";" _tuple ${_ds})
      list(GET _tuple 0 _name)
      list(GET _tuple 1 _vers)
      list(GET _tuple 2 _fnam)
      list(GET _tuple 3 _suffix)

      ExternalProject_Add(${_name} 
        URL ${_urlprefix}/${_fnam}.${_vers}.${_suffix}
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
endif()



