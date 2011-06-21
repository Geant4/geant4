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
# At present we provide the option, but it doesn't do anything if enabled
#
option(GEANT4_INSTALL_DATA "Install Geant4 Data Libraries (NOTIMPLEMENTED)" OFF)

# DO NOTHING YET...


