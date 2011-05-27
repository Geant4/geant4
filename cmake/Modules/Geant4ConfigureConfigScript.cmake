# - Script for configuring and installing geant4-config script
#
# The geant4-config script provides an sh based interface to provide
# information on the Geant4 installation, including installation prefix,
# version number, compiler and linker flags.
#
# The script is generated froma template file and then installed to the
# known bindir as an executable.
#
# $Id: Geant4ConfigureConfigScript.cmake,v 1.4 2010-12-13 17:31:59 bmorgan Exp $
# GEANT4 Tag $Name: not supported by cvs2svn $
#

#-----------------------------------------------------------------------------
# Only create script if we have a global library build...
#
if(NOT GEANT4_BUILD_GRANULAR_LIBS AND UNIX)
    # Setup variables needed for expansion in configuration file
    # - GDML
    if(GEANT4_USE_GDML)
        set(G4_BUILTWITH_GDML "yes")
    else()
        set(G4_BUILTWITH_GDML "no")
    endif()

    # - Qt
    if(GEANT4_USE_QT)
        set(G4_BUILTWITH_QT "yes")
    else()
        set(G4_BUILTWITH_QT "no")
    endif()

    # - RayTracerX
    if(GEANT4_USE_RAYTRACERX)
        set(G4_BUILTWITH_RAYTRACERX11 "yes")
        # We have to play with the X11 paths to get a clean set suitable for
        # inclusion
        set(_raw_x11_includes ${X11_INCLUDE_DIR})
        list(REMOVE_DUPLICATES _raw_x11_includes)
        set(G4_X11_INCLUDE_STATEMENT )
        foreach(_p ${_raw_x11_includes})
            set(G4_X11_INCLUDE_STATEMENT "-I${_p} ${G4_X11_INCLUDE_STATEMENT}")
        endforeach()
    else()
        set(G4_BUILTWITH_RAYTRACERX11 "no")
    endif()

    # Configure the script
    configure_file(${CMAKE_SOURCE_DIR}/cmake/Templates/geant4-config.in
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/geant4-config
        @ONLY)

    # Install it
    install(FILES ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/geant4-config
        DESTINATION ${GEANT4_BINDIR}
        PERMISSIONS
            OWNER_READ OWNER_WRITE OWNER_EXECUTE
            GROUP_READ GROUP_EXECUTE
            WORLD_READ WORLD_EXECUTE)
else()
    message(WARNING "geant4-config script will not be generated")
endif()

