# - Define useful macros for building and installing Geant4 library targets
#
# This file defines the following macros for Geant4 developers needing to
# add shared and static library targets.
#
# GEANT4_LIBRARY_TARGET        - define standard Geant4 library targets
#
# The macro will take the name of the library and its sources, defining
# static and shared targets depending on the value of BUILD_SHARED_LIBS and
# BUILD_STATIC_LIBS. Install targets are also created.

if(__GEANT4MACROLIBRARYTARGETS_ISLOADED)
  return()
endif()
set(__GEANT4MACROLIBRARYTARGETS_ISLOADED TRUE)

include(CMakeMacroParseArguments)

#----------------------------------------------------------------------------
# - GEANT4_LIBRARY_TARGET
# General build and install of a Geant4 library target
#
MACRO(GEANT4_LIBRARY_TARGET)
   CMAKE_PARSE_ARGUMENTS(G4LIBTARGET
        ""
        "NAME" "SOURCES;GEANT4_LINK_LIBRARIES;LINK_LIBRARIES"
        ${ARGN}
    )

    if(BUILD_SHARED_LIBS)
        # Add the shared library target and link its dependencies
        # WIN32 first
        if(WIN32)
            # We have to generate the def export file from an archive library.
            # If we're building Static libraries already, use that existing
            # target, otherwise, build a temporary uninstalled archive...
            if(BUILD_STATIC_LIBS)
                set(_archive ${G4LIBTARGET_NAME}-static)
            else()
                add_library(_${G4LIBTARGET_NAME}-archive STATIC EXCLUDE_FROM_ALL ${G4LIBTARGET_SOURCES})
                set(_archive _${G4LIBTARGET_NAME}-archive)
            endif()

            # - Create the .def file for this library
            # Note that we have to pass the actual full path to the library
            # to the command. CMake unfortunately won't generate this for us.
            # Note also that because we're likely to be on a platform with
            # multiconfig build tools. we use the CFG_INTDIR to locate the
            # archive we need...
            add_custom_command(OUTPUT _${G4LIBTARGET_NAME}.def
                COMMAND genwindef -o _${G4LIBTARGET_NAME}.def -l ${G4LIBTARGET_NAME} ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/${CMAKE_CFG_INTDIR}/${_archive}.lib
                DEPENDS ${_archive} genwindef)

            # - Now we can build the DLL
            # We create it from a dummy empty C++ file plus the def file.
            file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/_${G4LIBTARGET_NAME}.cpp
                "// empty _${G4LIBTARGET_NAME}.cpp\n")

            add_library(${G4LIBTARGET_NAME} SHARED _${G4LIBTARGET_NAME}.cpp
                _${G4LIBTARGET_NAME}.def)

            # - Link the DLL.
            # We link it to the archive, and the supplied libraries,
            # but then remove the archive from the LINK_INTERFACE.
            target_link_libraries(${G4LIBTARGET_NAME}
                ${_archive}
                ${G4LIBTARGET_GEANT4_LINK_LIBRARIES} 
                ${G4LIBTARGET_LINK_LIBRARIES})

            set_target_properties(${G4LIBTARGET_NAME}
                PROPERTIES LINK_INTERFACE_LIBRARIES "${G4LIBTARGET_GEANT4_LINK_LIBRARIES};${G4LIBTARGET_LINK_LIBRARIES}")

        else()
            # - We build a Shared library in the usual fashion...
            add_library(${G4LIBTARGET_NAME} SHARED ${G4LIBTARGET_SOURCES})
            target_link_libraries(${G4LIBTARGET_NAME}
                ${G4LIBTARGET_GEANT4_LINK_LIBRARIES} 
                ${G4LIBTARGET_LINK_LIBRARIES})
        endif()

        # This property is set to prevent concurrent builds of static and shared
        # libs removing each others files.
        set_target_properties(${G4LIBTARGET_NAME} 
            PROPERTIES CLEAN_DIRECT_OUTPUT 1)

        # Set the INSTALL_NAME_DIR of the library to its final installation
        # location (Only affects Mac OS X). This will only affect the library
        # when installed, BUT it does hard code this in. One should still be
        # able to bundle up the libraries later as CMake should build the
        # library with headerpad_max_install_names
        set_target_properties(${G4LIBTARGET_NAME}
            PROPERTIES INSTALL_NAME_DIR ${CMAKE_INSTALL_FULL_LIBDIR})

        # Install the library - note the use of RUNTIME, LIBRARY and ARCHIVE
        # this helps with later DLL builds.
        # Export to standard depends file for later install
        # NEEDS WORK TO REMOVE HARDCODED LIB/BIN DIR
        install(TARGETS ${G4LIBTARGET_NAME}
            EXPORT Geant4LibraryDepends
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

        # Append the library target to a global property so that build tree
        # export of library dependencies can pick up all targets
        set_property(GLOBAL APPEND
            PROPERTY GEANT4_EXPORTED_TARGETS ${G4LIBTARGET_NAME})
    endif()

    #
    # As above, but for static rather than shared library
    if(BUILD_STATIC_LIBS)
        # We have to distinguish the static from shared lib, so use -static in
        # name. Link its dependencies, and ensure we actually link to the
        # -static targets (We should strictly do this for the external
        # libraries as well if we want a pure static build).
        add_library(${G4LIBTARGET_NAME}-static STATIC ${G4LIBTARGET_SOURCES})

        set(G4LIBTARGET_GEANT4_LINK_LIBRARIES_STATIC )
        foreach(_tgt ${G4LIBTARGET_GEANT4_LINK_LIBRARIES})
            list(APPEND G4LIBTARGET_GEANT4_LINK_LIBRARIES_STATIC ${_tgt}-static)
        endforeach()

        target_link_libraries(${G4LIBTARGET_NAME}-static 
            ${G4LIBTARGET_GEANT4_LINK_LIBRARIES_STATIC}
            ${G4LIBTARGET_LINK_LIBRARIES})

        # But we can rename the output library to the correct name
        # On WIN32 we *retain* the -static postfix because otherwise
        # we'll conflict with the .lib from the DLL build...
        # We could also install differently...
        if(NOT WIN32)
            set_target_properties(${G4LIBTARGET_NAME}-static 
                PROPERTIES OUTPUT_NAME ${G4LIBTARGET_NAME})
        endif()
    
        set_target_properties(${G4LIBTARGET_NAME}-static
        	PROPERTIES CLEAN_DIRECT_OUTPUT 1)

        install(TARGETS ${G4LIBTARGET_NAME}-static
            EXPORT Geant4LibraryDepends
        	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
	        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

        set_property(GLOBAL APPEND
            PROPERTY GEANT4_EXPORTED_TARGETS ${G4LIBTARGET_NAME}-static)
    endif()
ENDMACRO()


#----------------------------------------------------------------------------
# - GEANT4_HEADER_MODULE_TARGET
# Build and install for a header only Geant4 module.
#
MACRO(GEANT4_HEADER_MODULE_TARGET)
    CMAKE_PARSE_ARGUMENTS(G4HEADERMOD
        ""
        "COMPONENT"
        ""
        ${ARGN}
    )

    # Only has one component, and we just have to pick out the headers
    include(${G4HEADERMOD_COMPONENT})

    # Header install?
    install(FILES ${${G4MODULENAME}_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
        COMPONENT Development)

    # Store the include path of the component so that the build tree
    # config file can pick up all needed header paths
    set_property(GLOBAL APPEND
        PROPERTY GEANT4_BUILDTREE_INCLUDE_DIRS ${${G4MODULENAME}_INCDIR})
ENDMACRO()



#----------------------------------------------------------------------------
# - GEANT4_GRANULAR_LIBRARY_TARGET
# Build and install for a Geant4 module (granular) library
#
MACRO(GEANT4_GRANULAR_LIBRARY_TARGET)
    CMAKE_PARSE_ARGUMENTS(G4GRANLIB
        ""
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
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
        COMPONENT Development)

    # Store the include path of the component so that the build tree
    # config file can pick up all needed header paths
    set_property(GLOBAL APPEND
        PROPERTY GEANT4_BUILDTREE_INCLUDE_DIRS ${${G4MODULENAME}_INCDIR})
ENDMACRO()



#----------------------------------------------------------------------------
# - GEANT4_GLOBAL_LIBRARY_TARGET
# Build and install of a Geant4 category (global) library
#
MACRO(GEANT4_GLOBAL_LIBRARY_TARGET)
    CMAKE_PARSE_ARGUMENTS(G4GLOBLIB
        ""
        "NAME"
        "COMPONENTS"
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

        list(APPEND ${G4GLOBLIB_NAME}_BUILDTREE_INCLUDES ${${G4MODULENAME}_INCDIR})
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
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
        COMPONENT Development)

    # Store the include path of the component so that the build tree
    # config file can pick up all needed header paths
    set_property(GLOBAL APPEND
        PROPERTY GEANT4_BUILDTREE_INCLUDE_DIRS ${${G4GLOBLIB_NAME}_BUILDTREE_INCLUDES})

ENDMACRO()



