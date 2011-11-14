# - Script for configuring and installing geant4-config script
#
# The geant4-config script provides an sh based interface to provide
# information on the Geant4 installation, including installation prefix,
# version number, compiler and linker flags.
#
# The script is generated from a template file and then installed to the
# known bindir as an executable.
#
# Paths are always hardcoded in the build tree version as this is never
# intended to be relocatable.
# The Install Tree script uses self-location based on that in 
# {root,clehep}-config is the install itself is relocatable, otherwise
# absolute paths are encoded.
#
# $Id: Geant4ConfigureConfigScript.cmake,v 1.4 2010-12-13 17:31:59 bmorgan Exp $
# GEANT4 Tag $Name: not supported by cvs2svn $
#

#-----------------------------------------------------------------------------
# Only create script if we have a global library build...
#
if(NOT GEANT4_BUILD_GRANULAR_LIBS AND UNIX)
  # Setup variables needed for expansion in configuration file
  # - CLHEP
  if(GEANT4_USE_SYSTEM_CLHEP)
    set(G4_BUILTWITH_CLHEP "no")
  else()
    set(G4_BUILTWITH_CLHEP "yes")
  endif()

  # - EXPAT
  if(GEANT4_USE_SYSTEM_EXPAT)
    set(G4_BUILTWITH_EXPAT "no")
  else()
    set(G4_BUILTWITH_EXPAT "yes")
  endif()

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

  # - Motif
  if(GEANT4_USE_XM)
    set(G4_BUILTWITH_MOTIF "yes")
    set(G4_CONFIG_NEEDS_X11 TRUE)
  else()
    set(G4_BUILTWITH_MOTIF "no")
  endif()

  # - RayTracerX
  if(GEANT4_USE_RAYTRACER_X11)
    set(G4_BUILTWITH_RAYTRACERX11 "yes")
    set(G4_CONFIG_NEEDS_X11 TRUE)
  else()
    set(G4_BUILTWITH_RAYTRACERX11 "no")
  endif()

  # - OpenGL X11
  if(GEANT4_USE_OPENGL_X11)
    set(G4_BUILTWITH_OPENGLX11 "yes")
    set(G4_CONFIG_NEEDS_X11 TRUE)
  else()
    set(G4_BUILTWITH_OPENGLX11 "no")
  endif()

  # - OpenInventor
  if(GEANT4_USE_INVENTOR)
    set(G4_BUILTWITH_INVENTOR "yes")
  else()
    set(G4_BUILTWITH_INVENTOR "no")
  endif()

  # If we have a module that uses X11, We have to play with the X11 paths to 
  # get a clean set suitable for inclusion
  if(G4_CONFIG_NEEDS_X11)
    set(_raw_x11_includes ${X11_INCLUDE_DIR})
    list(REMOVE_DUPLICATES _raw_x11_includes)
    set(G4_X11_INCLUDE_STATEMENT )
    foreach(_p ${_raw_x11_includes})
      set(G4_X11_INCLUDE_STATEMENT "-I${_p} ${G4_X11_INCLUDE_STATEMENT}")
    endforeach()
  endif()

  # Configure the script
  # - BUILD TREE
  # Ouch, the include path will be LONG, but at least we always have absolute
  # paths...
  set(GEANT4_CONFIG_SELF_LOCATION "# BUILD TREE IS NON-RELOCATABLE")
  set(GEANT4_CONFIG_INSTALL_PREFIX "${PROJECT_BINARY_DIR}")
  set(GEANT4_CONFIG_INSTALL_EXECPREFIX \"\")
  set(GEANT4_CONFIG_LIBDIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

  get_property(__geant4_buildtree_include_dirs GLOBAL PROPERTY
    GEANT4_BUILDTREE_INCLUDE_DIRS)

  foreach(_dir ${__geant4_buildtree_include_dirs})
    set(GEANT4_CONFIG_INCLUDE_DIRS "${GEANT4_CONFIG_INCLUDE_DIRS} \\
    ${_dir}")
  endforeach()

  # Configure the build tree script
  # If we're on CMake 2.8 and above, we try to use file(COPY) to create an
  # executable script
  # Not sure if version check is o.k., but I'll be shocked if we ever see
  # a CMake 2.7 in the wild...
  if(${CMAKE_VERSION} VERSION_GREATER 2.7)
    configure_file(${CMAKE_SOURCE_DIR}/cmake/Templates/geant4-config.in
      ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/geant4-config
      @ONLY)

    file(COPY 
      ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/geant4-config
      DESTINATION ${PROJECT_BINARY_DIR}
      FILE_PERMISSIONS
      OWNER_READ OWNER_WRITE OWNER_EXECUTE
      GROUP_READ GROUP_EXECUTE
      WORLD_READ WORLD_EXECUTE)
  else()
    # Changing permissions is awkward, so just configure and document
    # that you have to do 'sh geant4-config' in this case.
    configure_file(${CMAKE_SOURCE_DIR}/cmake/Templates/geant4-config.in
      ${PROJECT_BINARY_DIR}/geant4-config
      @ONLY
      )
  endif()

  # - Install Tree
  # Much easier :-)
  # Non-Relocatable case...
  if(CMAKE_INSTALL_IS_NONRELOCATABLE)
    # Hardcoded paths
    set(GEANT4_CONFIG_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")
    set(GEANT4_CONFIG_INSTALL_EXECPREFIX \"\")
    set(GEANT4_CONFIG_LIBDIR "${CMAKE_INSTALL_FULL_LIBDIR}")
    set(GEANT4_CONFIG_INCLUDE_DIRS "${CMAKE_INSTALL_FULL_INCLUDEDIR}/Geant4")
  else()
    # Calculate base of self contained install based on relative path from 
    # CMAKE_INSTALL_FULL_BINDIR to CMAKE_INSTALL_PREFIX.
    file(RELATIVE_PATH _bin_to_prefix ${CMAKE_INSTALL_FULL_BINDIR} ${CMAKE_INSTALL_PREFIX})
    # Strip any trailing path separators just for neatness.
    string(REGEX REPLACE "[/\\]$" "" _bin_to_prefix "${_bin_to_prefix}")

    set(GEANT4_CONFIG_INSTALL_PREFIX "$scriptloc/${_bin_to_prefix}")
    set(GEANT4_CONFIG_INSTALL_EXECPREFIX \"\")
    set(GEANT4_CONFIG_LIBDIR "\${prefix}/${CMAKE_INSTALL_LIBDIR}")
    set(GEANT4_CONFIG_INCLUDE_DIRS "\${prefix}/${CMAKE_INSTALL_INCLUDEDIR}/Geant4")
  endif() 

  # Configure the install tree script
  configure_file(${CMAKE_SOURCE_DIR}/cmake/Templates/geant4-config.in
    ${PROJECT_BINARY_DIR}/InstallTreeFiles/geant4-config
    @ONLY
    )

  # Install it
  install(FILES ${PROJECT_BINARY_DIR}/InstallTreeFiles/geant4-config
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    PERMISSIONS
      OWNER_READ OWNER_WRITE OWNER_EXECUTE
      GROUP_READ GROUP_EXECUTE
      WORLD_READ WORLD_EXECUTE
    COMPONENT Development
    )
endif()

