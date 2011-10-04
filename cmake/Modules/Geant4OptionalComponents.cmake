# - Setup core required and optional components of Geant4
#
# Here we provide options to enable and configure required and optional
# components, which may require third party libraries.
#
# We don't configure User Interface options here because these require
# a higher degree of configuration so to keep things neat these have their
# own Module.
#
# Options configured here:
#  CLHEP  - Control use of internal G4clhep, or locate external CLHEP
#  EXPAT  - Control use of internal G4expat, or locate external EXPAT.
#  ZLIB   - Control use of internal G4zlib, or locate external ZLIB
#           (NOTIMPLEMENTEDYET - always uses internal zlib)
#  GDML   - Requires external XercesC
#  G3TOG4 - UNIX only, requires fortran compiler


#----------------------------------------------------------------------------
# Find required CLHEP package
# We prefer to use our internal CLHEP library, but provide an option
# to allow an external CLHEP to be used should the user require this. 
# We also allow that it can be automatically enabled by providing
# the CLHEP_ROOT_DIR option (which FindCLHEP will recognize)
#
# KNOWNISSUE : For internal CLHEP, how to deal with static and shared?
if(CLHEP_ROOT_DIR)
    set(_default_use_system_clhep ON)
else()
    set(_default_use_system_clhep OFF)
endif()

option(GEANT4_USE_SYSTEM_CLHEP 
    "Use external CLHEP library" 
    ${_default_use_system_clhep}
)

if(GEANT4_USE_SYSTEM_CLHEP)
    # We keep this as required, because if the user chooses to use a
    # system option we assume that we absolutely, positively require this.
    find_package(CLHEP 2.1.0.1 REQUIRED)
    set(GEANT4_USE_SYSTEM_CLHEP TRUE)
else()
    set(CLHEP_FOUND TRUE)
    set(CLHEP_INCLUDE_DIRS
        "${PROJECT_SOURCE_DIR}/source/externals/clhep/include")
    if(BUILD_SHARED_LIBS)
        set(CLHEP_LIBRARIES G4clhep)
    else()
        set(CLHEP_LIBRARIES G4clhep-static)
    endif()
endif()

GEANT4_ADD_FEATURE(GEANT4_USE_SYSTEM_CLHEP "Using external install of CLHEP")


#------------------------------------------------------------------------------
# Find required EXPAT package
# We always use the internal G4expat on WIN32, on all other systems we prefer
# system install of expat, unless we can't find it. In that case, we fall 
# back to the internal G4expat
# If we use the internal G4expat, fix the EXPAT_XXX variables to point to the
# internal location of expat headers and library target so Geant4 clients of
# expat have a consistent interface.
#
# KNOWNISSUE : For internal expat, how to deal with static and shared?
if(WIN32)
    set(EXPAT_FOUND TRUE)
    set(EXPAT_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/source/externals/expat/include)
    set(EXPAT_LIBRARIES G4expat)
else()
    find_package(EXPAT)
    if(EXPAT_FOUND)
       set(GEANT4_USE_SYSTEM_EXPAT TRUE)
    else()
        set(EXPAT_FOUND TRUE)
        set(EXPAT_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/source/externals/expat/include)
        if(BUILD_SHARED_LIBS)
            set(EXPAT_LIBRARIES G4expat)
        else()
            set(EXPAT_LIBRARIES G4expat-static)
        endif()
    endif()
endif()

GEANT4_ADD_FEATURE(GEANT4_USE_SYSTEM_EXPAT "Using system install of EXPAT")


#------------------------------------------------------------------------------
# Find required ZLIB package
# For now, we always use the internal Geant4 zlib...
#
#option(GEANT4_USE_SYSTEM_ZLIB "Use the system's zlib library" OFF)
if(GEANT4_USE_SYSTEM_ZLIB)
    # This needs more work - use ITK's way of doing it as an example.
    find_package(ZLIB)
endif(GEANT4_USE_SYSTEM_ZLIB)

GEANT4_ADD_FEATURE(GEANT4_USE_SYSTEM_ZLIB "Use system zlib library")



#----------------------------------------------------------------------------
# Optional Support for GDML - requires Xerces-C package
#
if(XERCESC_ROOT_DIR)
    set(_default_use_gdml ON)
else()
    set(_default_use_gdml OFF)
endif()

option(GEANT4_USE_GDML
    "Build Geant4 with GDML support" 
    ${_default_use_gdml}
)

if(GEANT4_USE_GDML)
    find_package(XercesC REQUIRED)
endif(GEANT4_USE_GDML)

GEANT4_ADD_FEATURE(GEANT4_USE_GDML "Build Geant4 with GDML support")


#----------------------------------------------------------------------------
# Optional support for G3TOG4 convertion interface.
# We do not build the rztog4 application.
# -- OLDER NOTES --
# The G3toG4 *library* should always be built, but the rztog4 application 
# requires a Fortran compiler AND CERNLIB, so is optional.
# Only on *NIX because converter requires CERNLIB, and Windows support for
# this is an unknown quantity at present (can always change later).
# -- OLDER NOTES --
#
if(UNIX)
    option(GEANT4_USE_G3TOG4 "Build the rztog4 application for converting Geant3 geometries to Geant4" OFF)
    #if(GEANT4_USE_G3TOG4)
    #enable_language(Fortran)
    #find_package(CERNLIB REQUIRED)
    #endif(GEANT4_USE_G3TOG4)

    GEANT4_ADD_FEATURE(GEANT4_USE_G3TOG4 "Build the rztog4 application")
endif()









