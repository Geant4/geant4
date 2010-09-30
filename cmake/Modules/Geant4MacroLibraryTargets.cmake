# - Define useful macros for building and installing library targets
#
# This file defines the following macros for Geant4 developers needing to
# add shared and static library targets.
#
# GEANT4_LIBRARY_TARGET        - define standard Geant4 library targets
#
# The macro will take the name of the library and its sources, defining
# static and shared targets depending on the value of BUILD_SHARED_LIBS and
# BUILD_STATIC_LIBS. Install targets are also created.

include(CMakeMacroParseArguments)

MACRO(GEANT4_LIBRARY_TARGET)
    PARSE_ARGUMENTS(G4LIBTARGET
        "NAME;SOURCES;GEANT4_LINK_LIBRARIES;LINK_LIBRARIES"
        ""
        ${ARGN}
    )

    if(BUILD_SHARED_LIBS)
        # Add the shared library target and link its dependencies
        add_library(${G4LIBTARGET_NAME} SHARED ${G4LIBTARGET_SOURCES})
        target_link_libraries(${G4LIBTARGET_NAME}
            ${G4LIBTARGET_GEANT4_LINK_LIBRARIES} 
            ${G4LIBTARGET_LINK_LIBRARIES})

        # This property is set to prevent concurrent builds of static and shared
        # libs removing each others files.
        set_target_properties(${G4LIBTARGET_NAME} 
            PROPERTIES CLEAN_DIRECT_OUTPUT 1)

        # Install the library - note the use of RUNTIME, LIBRARY and ARCHIVE
        # this helps with later DLL builds
        # NEEDS WORK TO REMOVE HARDCODED LIB/BIN DIR
        install(TARGETS ${G4LIBTARGET_NAME}
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)
    endif()

    #
    # As above, but for static rather than shared library
    if(BUILD_STATIC_LIBS)
        # We have to distinguish the static from shared lib, so use -static in
        # name. Link its dependencies, but CMake should only use this for
        # tracking.
        add_library(${G4LIBTARGET_NAME}-static STATIC ${G4LIBTARGET_SOURCES})
        target_link_libraries(${G4LIBTARGET_NAME}-static 
            ${G4LIBTARGET_GEANT4_LINK_LIBRARIES}
            ${G4LIBTARGET_LINK_LIBRARIES})

        #But we can rename the output library to the correct name
        set_target_properties(${G4LIBTARGET_NAME}-static 
            PROPERTIES OUTPUT_NAME ${G4LIBTARGET_NAME})
    
        set_target_properties(${G4LIBTARGET_NAME}-static
        	PROPERTIES CLEAN_DIRECT_OUTPUT 1)

        install(TARGETS ${G4LIBTARGET_NAME}-static
        	RUNTIME DESTINATION bin
	        LIBRARY DESTINATION lib
        	ARCHIVE DESTINATION lib)
    endif()
ENDMACRO()


MACRO(GEANT4_HEADER_MODULE_TARGET)
    PARSE_ARGUMENTS(G4HEADERMOD
        "COMPONENT"
        ""
        ${ARGN}
    )

    # Only has one component, and we just have to pick out the headers
    include(${G4HEADERMOD_COMPONENT})

    # Header install?
    install(FILES ${${G4MODULENAME}_HEADERS}
        DESTINATION include/${CMAKE_PROJECT_NAME}
        COMPONENT headers)
ENDMACRO()


MACRO(GEANT4_GRANULAR_LIBRARY_TARGET)
    PARSE_ARGUMENTS(G4GRANLIB
        "COMPONENT"
        ""
        ${ARGN}
    )

    # Granular lib only has one component, but we must pick out
    # the granular dependencies
    include(${G4GRANLIB_COMPONENT})

    # Add the library target, using variables set by the inclusion of
    # the component file
    GEANT4_LIBRARY_TARGET(NAME ${G4MODULENAME}
        SOURCES ${${G4MODULENAME}_SOURCES} ${${G4MODULENAME}_HEADERS}
        GEANT4_LINK_LIBRARIES ${${G4MODULENAME}_GRANULAR_DEPENDENCIES}
        LINK_LIBRARIES ${${G4MODULENAME}_LINK_LIBRARIES})

    # Header install?
    install(FILES ${${G4MODULENAME}_HEADERS}
        DESTINATION include/${CMAKE_PROJECT_NAME}
        COMPONENT headers)
ENDMACRO()



MACRO(GEANT4_GLOBAL_LIBRARY_TARGET)
    PARSE_ARGUMENTS(G4GLOBLIB
        "NAME;COMPONENTS"
        ""
        ${ARGN}
    )

    # We loop over the component sources one at a time,
    # appending properties as we go.
    foreach(_comp ${G4GLOBLIB_COMPONENTS})
        include(${_comp})
        # In case we have a global lib with one component, ensure name gets set
        if(NOT G4GLOBLIB_NAME)
            set(G4GLOBLIB_NAME ${G4MODULENAME})
        endif()

        list(APPEND ${G4GLOBLIB_NAME}_GLOBAL_SOURCES ${${G4MODULENAME}_SOURCES})
        list(APPEND ${G4GLOBLIB_NAME}_GLOBAL_HEADERS ${${G4MODULENAME}_HEADERS})

        list(APPEND ${G4GLOBLIB_NAME}_GLOBAL_DEPENDENCIES
            ${${G4MODULENAME}_GLOBAL_DEPENDENCIES})

        list(APPEND ${G4GLOBLIB_NAME}_LINK_LIBRARIES
            ${${G4MODULENAME}_LINK_LIBRARIES})
    endforeach()

    # Filter out duplicates in GLOBAL_DEPENDENCIES and LINK_LIBRARIES
    if(${G4GLOBLIB_NAME}_GLOBAL_DEPENDENCIES)
        list(REMOVE_DUPLICATES ${G4GLOBLIB_NAME}_GLOBAL_DEPENDENCIES)
    endif()
    if(${G4GLOBLIB_NAME}_LINK_LIBRARIES)
        list(REMOVE_DUPLICATES ${G4GLOBLIB_NAME}_LINK_LIBRARIES)
    endif()

    # Now add the library target
    GEANT4_LIBRARY_TARGET(NAME ${G4GLOBLIB_NAME} 
        SOURCES 
            ${${G4GLOBLIB_NAME}_GLOBAL_SOURCES} 
            ${${G4GLOBLIB_NAME}_GLOBAL_HEADERS}
        GEANT4_LINK_LIBRARIES 
            ${${G4GLOBLIB_NAME}_GLOBAL_DEPENDENCIES}
        LINK_LIBRARIES 
            ${${G4GLOBLIB_NAME}_LINK_LIBRARIES})

    # Header install?
    install(FILES ${${G4GLOBLIB_NAME}_GLOBAL_HEADERS}
        DESTINATION include/${CMAKE_PROJECT_NAME}
        COMPONENT headers)
ENDMACRO()



