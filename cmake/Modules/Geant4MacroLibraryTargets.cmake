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
#
# A custom compile definition "GEANT4_DEVELOPER_<CONFIG>" is set on
# each target using the target property COMPILE_DEFINITIONS_<CONFIG>
# target property

if(__GEANT4MACROLIBRARYTARGETS_ISLOADED)
  return()
endif()
set(__GEANT4MACROLIBRARYTARGETS_ISLOADED TRUE)

include(CMakeMacroParseArguments)

#-----------------------------------------------------------------------
# function geant4_compile_definitions_config(<target>)
#          Set a custom compile definition for a Geant4 target on a
#          per configuration basis:
#            For mode <CONFIG>, define GEANT4_DEVELOPER_<CONFIG>
#
function(geant4_compile_definitions_config _target)
  if(NOT TARGET ${_target})
    message(FATAL_ERROR "geant4_compile_definitions_config passed target '${_target}' which is not a valid CMake target")
  endif()

  # CMake 3 prefers $<CONFIG> in generator expressions, but 2.8.12 only
  # supports $<CONFIGURATION>
  set(__default_config_generator "\$<CONFIG>")
  if(CMAKE_MAJOR_VERSION LESS 3)
    set(__default_config_generator "\$<CONFIGURATION>")
  endif()

  set_property(TARGET ${_target}
    APPEND PROPERTY COMPILE_DEFINITIONS
      GEANT4_DEVELOPER_${__default_config_generator}
    )
endfunction()


#-----------------------------------------------------------------------
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
      # This is a temporary separate from a real archive library, and
      # even though it's static, we need to mark that it will have
      # DLL symbols via the G4LIB_BUILD_DLL macro
      add_library(_${G4LIBTARGET_NAME}-archive STATIC EXCLUDE_FROM_ALL ${G4LIBTARGET_SOURCES})
      set(_archive _${G4LIBTARGET_NAME}-archive)
      target_compile_features(${_archive} PUBLIC ${GEANT4_TARGET_COMPILE_FEATURES})
      target_compile_definitions(${_archive} PUBLIC -DG4LIB_BUILD_DLL)

      # - Add the config specific compile definitions
      geant4_compile_definitions_config(${_archive})

      # - Create the .def file for this library
      # Use generator expressions to get path to per-mode lib and
      # older CMAKE_CFG_INTDIR variable to set name of per-mode def
      # file (Needed as generator expressions cannot be used in argument
      # to OUTPUT...
      add_custom_command(OUTPUT _${G4LIBTARGET_NAME}-${CMAKE_CFG_INTDIR}.def
        COMMAND genwindef -o _${G4LIBTARGET_NAME}-${CMAKE_CFG_INTDIR}.def -l ${G4LIBTARGET_NAME} $<TARGET_FILE:${_archive}>
        DEPENDS ${_archive} genwindef)

      # - Now we can build the DLL
      # We create it from a dummy empty C++ file plus the def file.
      # Also set the public compile definition on it so that clients
      # will set correct macro automatically.
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/_${G4LIBTARGET_NAME}.cpp
        "// empty _${G4LIBTARGET_NAME}.cpp\n")

      add_library(${G4LIBTARGET_NAME} SHARED _${G4LIBTARGET_NAME}.cpp
        _${G4LIBTARGET_NAME}-${CMAKE_CFG_INTDIR}.def)
      target_compile_definitions(${G4LIBTARGET_NAME} PUBLIC -DG4LIB_BUILD_DLL)
      target_compile_features(${G4LIBTARGET_NAME} PUBLIC ${GEANT4_TARGET_COMPILE_FEATURES})

      # - Link the DLL.
      # We link it to the archive, and the supplied libraries,
      # but then remove the archive from the LINK_INTERFACE.
      target_link_libraries(${G4LIBTARGET_NAME}
        ${_archive}
        ${G4LIBTARGET_GEANT4_LINK_LIBRARIES}
        ${G4LIBTARGET_LINK_LIBRARIES})

      set_target_properties(${G4LIBTARGET_NAME}
        PROPERTIES INTERFACE_LINK_LIBRARIES "${G4LIBTARGET_GEANT4_LINK_LIBRARIES};${G4LIBTARGET_LINK_LIBRARIES}")

    else()
      # - We build a Shared library in the usual fashion...
      add_library(${G4LIBTARGET_NAME} SHARED ${G4LIBTARGET_SOURCES})
      geant4_compile_definitions_config(${G4LIBTARGET_NAME})
      target_compile_features(${G4LIBTARGET_NAME} PUBLIC ${GEANT4_TARGET_COMPILE_FEATURES})
      target_link_libraries(${G4LIBTARGET_NAME}
        ${G4LIBTARGET_GEANT4_LINK_LIBRARIES}
        ${G4LIBTARGET_LINK_LIBRARIES})
    endif()

    # This property is set to prevent concurrent builds of static and
    # shared libs removing each others files.
    set_target_properties(${G4LIBTARGET_NAME}
      PROPERTIES CLEAN_DIRECT_OUTPUT 1)

    # Always use '@rpath' in install names of libraries. This is the
    # most flexible mechanism for
    set_target_properties(${G4LIBTARGET_NAME}
      PROPERTIES MACOSX_RPATH 1
      )

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
    geant4_compile_definitions_config(${G4LIBTARGET_NAME}-static)
    target_compile_features(${G4LIBTARGET_NAME}-static PUBLIC ${GEANT4_TARGET_COMPILE_FEATURES})

    set(G4LIBTARGET_GEANT4_LINK_LIBRARIES_STATIC )
    foreach(_tgt ${G4LIBTARGET_GEANT4_LINK_LIBRARIES})
      list(APPEND G4LIBTARGET_GEANT4_LINK_LIBRARIES_STATIC ${_tgt}-static)
    endforeach()

    # If we are building both types of library and builtin clhep etc,
    # we want to link shared->shared and static->static.
    # Because externals like clhep appear in G4LIBTARGET_LINK_LIBRARIES,
    # filter this list to replace shared builtins with their static variant
    string(REGEX REPLACE
      "(G4clhep|G4expat|G4zlib|G4geomUSolids)(;|$)" "\\1-static\\2"
      G4LIBTARGET_LINK_LIBRARIES_STATIC
      "${G4LIBTARGET_LINK_LIBRARIES}"
      )

    target_link_libraries(${G4LIBTARGET_NAME}-static
      ${G4LIBTARGET_GEANT4_LINK_LIBRARIES_STATIC}
      ${G4LIBTARGET_LINK_LIBRARIES_STATIC})

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

#-----------------------------------------------------------------------
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

#-----------------------------------------------------------------------
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

#-----------------------------------------------------------------------
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

  # Filter out duplicates/self in GLOBAL_DEPENDENCIES and LINK_LIBRARIES
  if(${G4GLOBLIB_NAME}_GLOBAL_DEPENDENCIES)
    list(REMOVE_DUPLICATES ${G4GLOBLIB_NAME}_GLOBAL_DEPENDENCIES)
    list(REMOVE_ITEM ${G4GLOBLIB_NAME}_GLOBAL_DEPENDENCIES ${G4GLOBLIB_NAME})
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



