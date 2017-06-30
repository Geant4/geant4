#.rst:
# G4ConfigurePkgConfigHelpers
# ---------------------------
#
# This module configures and installs pkg-config scripts and the
# geant4-config program to help clients compile and link against
# the Geant4 libraries.
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
# {root,clhep}-config is the install itself is relocatable, otherwise
# absolute paths are encoded.
#
#

#-----------------------------------------------------------------
# License and Disclaimer
#
# The  Geant4 software  is  copyright of the Copyright Holders  of
# the Geant4 Collaboration.  It is provided  under  the terms  and
# conditions of the Geant4 Software License,  included in the file
# LICENSE and available at  http://cern.ch/geant4/license .  These
# include a list of copyright holders.
#
# Neither the authors of this software system, nor their employing
# institutes,nor the agencies providing financial support for this
# work  make  any representation or  warranty, express or implied,
# regarding  this  software system or assume any liability for its
# use.  Please see the license in the file  LICENSE  and URL above
# for the full disclaimer and the limitation of liability.
#
# This  code  implementation is the result of  the  scientific and
# technical work of the GEANT4 collaboration.
# By using,  copying,  modifying or  distributing the software (or
# any work based  on the software)  you  agree  to acknowledge its
# use  in  resulting  scientific  publications,  and indicate your
# acceptance of all terms of the Geant4 Software license.
#
#-----------------------------------------------------------------


#-----------------------------------------------------------------------
# function get_system_include_dirs
#          return list of directories our C++ compiler searches
#          by default.
#
#          The idea comes from CMake's inbuilt technique to do this
#          for the Eclipse and CodeBlocks generators, but we implement
#          our own function because the CMake functionality is internal
#          so we can't rely on it.
function(get_system_include_dirs _dirs)
  # Only for GCC, Clang and Intel
  if("${CMAKE_CXX_COMPILER_ID}" MATCHES GNU OR "${CMAKE_CXX_COMPILER_ID}" MATCHES Clang OR "${CMAKE_CXX_COMPILER_ID}" MATCHES Intel)
    # Proceed
    file(WRITE "${CMAKE_BINARY_DIR}/CMakeFiles/g4dummy" "\n")

    # Save locale, them to "C" english locale so we can parse in English
    set(_orig_lc_all      $ENV{LC_ALL})
    set(_orig_lc_messages $ENV{LC_MESSAGES})
    set(_orig_lang        $ENV{LANG})

    set(ENV{LC_ALL}      C)
    set(ENV{LC_MESSAGES} C)
    set(ENV{LANG}        C)

    execute_process(COMMAND ${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ARG1} -v -E -x c++ -dD g4dummy
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles
      ERROR_VARIABLE _cxxOutput
      OUTPUT_VARIABLE _cxxStdout
      )

    file(REMOVE "${CMAKE_BINARY_DIR}/CMakeFiles/g4dummy")

    # Parse and extract search dirs
    set(_resultIncludeDirs )
    if( "${_cxxOutput}" MATCHES "> search starts here[^\n]+\n *(.+ *\n) *End of (search) list" )
      string(REGEX MATCHALL "[^\n]+\n" _includeLines "${CMAKE_MATCH_1}")
      foreach(nextLine ${_includeLines})
        string(REGEX REPLACE "\\(framework directory\\)" "" nextLineNoFramework "${nextLine}")
        string(STRIP "${nextLineNoFramework}" _includePath)
        list(APPEND _resultIncludeDirs "${_includePath}")
      endforeach()
    endif()

    # Restore original locale
    set(ENV{LC_ALL}      ${_orig_lc_all})
    set(ENV{LC_MESSAGES} ${_orig_lc_messages})
    set(ENV{LANG}        ${_orig_lang})

    set(${_dirs} ${_resultIncludeDirs} PARENT_SCOPE)
  else()
    set(${_dirs} "" PARENT_SCOPE)
  endif()
endfunction()

#-----------------------------------------------------------------------
# Only create script if we have a global library build...
#
if(NOT GEANT4_BUILD_GRANULAR_LIBS AND UNIX)
  # Get implicit search paths
  get_system_include_dirs(_cxx_compiler_dirs)

  # Setup variables needed for expansion in configuration file
  # - Static libs
  if(BUILD_STATIC_LIBS)
    set(G4_BUILTWITH_STATICLIBS "yes")
  else()
    set(G4_BUILTWITH_STATICLIBS "no")
  endif()

  # - Multithreading
  if(GEANT4_BUILD_MULTITHREADED)
    set(G4_BUILTWITH_MULTITHREADING "yes")
  else()
    set(G4_BUILTWITH_MULTITHREADING "no")
  endif()

  # - CLHEP
  if(GEANT4_USE_SYSTEM_CLHEP)
    set(G4_BUILTWITH_CLHEP "no")
    #inc path
    get_filename_component(G4_SYSTEM_CLHEP_INCLUDE_DIR "${CLHEP_INCLUDE_DIR}" ABSOLUTE)

    #libpath
    list(GET CLHEP_LIBRARIES 0 _zeroth_clhep_lib)
    get_target_property(_system_clhep_libdir "${_zeroth_clhep_lib}" LOCATION)
    get_filename_component(_system_clhep_libdir "${_system_clhep_libdir}" REALPATH)
    get_filename_component(_system_clhep_libdir "${_system_clhep_libdir}" DIRECTORY)
    set(G4_SYSTEM_CLHEP_LIBRARIES "-L${_system_clhep_libdir}")

    foreach(_clhep_lib ${CLHEP_LIBRARIES})
      get_target_property(_curlib "${_clhep_lib}" LOCATION)
      get_filename_component(_curlib "${_curlib}" NAME)
      string(REGEX REPLACE "^lib(.*)\\.(so|a|dylib|lib|dll)$" "\\1" _curlib "${_curlib}")
      set(G4_SYSTEM_CLHEP_LIBRARIES "${G4_SYSTEM_CLHEP_LIBRARIES} -l${_curlib}")
    endforeach()
  else()
    set(G4_BUILTWITH_CLHEP "yes")
  endif()

  # - EXPAT
  if(GEANT4_USE_SYSTEM_EXPAT)
    set(G4_BUILTWITH_EXPAT "no")
  else()
    set(G4_BUILTWITH_EXPAT "yes")
  endif()

  # - ZLIB
  if(GEANT4_USE_SYSTEM_ZLIB)
    set(G4_BUILTWITH_ZLIB "no")
  else()
    set(G4_BUILTWITH_ZLIB "yes")
  endif()

  # - GDML
  if(GEANT4_USE_GDML)
    set(G4_BUILTWITH_GDML "yes")
    set(G4_XERCESC_INCLUDE_DIRS ${XERCESC_INCLUDE_DIRS})
    list(REMOVE_DUPLICATES G4_XERCESC_INCLUDE_DIRS)
    list(REMOVE_ITEM G4_XERCESC_INCLUDE_DIRS ${_cxx_compiler_dirs})

    set(G4_XERCESC_CFLAGS )
    foreach(_dir ${G4_XERCESC_INCLUDE_DIRS})
      set(G4_XERCESC_CFLAGS "${G4_XERCESC_CFLAGS} -I${_dir}")
    endforeach()
  else()
    set(G4_BUILTWITH_GDML "no")
  endif()

  # - G3ToG4
  if(GEANT4_USE_G3TOG4)
    set(G4_BUILTWITH_G3TOG4 "yes")
  else()
    set(G4_BUILTWITH_G3TOG4 "no")
  endif()

  # - USolids
  if(GEANT4_USE_USOLIDS OR GEANT4_USE_PARTIAL_USOLIDS)
    set(G4_BUILTWITH_USOLIDS "yes")
    set(G4_USOLIDS_INCLUDE_DIRS "${USOLIDS_INCLUDE_DIRS} ${VECGEOM_EXTERNAL_INCLUDES}")
    list(REMOVE_DUPLICATES G4_USOLIDS_INCLUDE_DIRS)
    list(REMOVE_ITEM G4_USOLIDS_INCLUDE_DIRS ${_cxx_compiler_dirs})

    string(REPLACE ";" " " G4_USOLIDS_CFLAGS "${GEANT4_USOLIDS_COMPILE_DEFINITIONS}")
    foreach(_dir ${G4_USOLIDS_INCLUDE_DIRS})
      set(G4_USOLIDS_CFLAGS "${G4_USOLIDS_CFLAGS} -I${_dir}")
    endforeach()
  else()
    set(G4_BUILTWITH_USOLIDS "no")
  endif()

  # - Freetype
  if(GEANT4_USE_FREETYPE)
    set(G4_BUILTWITH_FREETYPE "yes")
  else()
    set(G4_BUILTWITH_FREETYPE "no")
  endif()

  # - Qt
  if(GEANT4_USE_QT)
    set(G4_BUILTWITH_QT "yes")
    if(QT4_FOUND)
      set(G4_QT_INCLUDE_DIRS ${QT_QTCORE_INCLUDE_DIR} ${QT_QTGUI_INCLUDE_DIR} ${QT_QTOPENGL_INCLUDE_DIR})
    else()
      set(G4_QT_INCLUDE_DIRS ${Qt5Core_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS} ${Qt5OpenGL_INCLUDE_DIRS} ${Qt5PrintSupport_INCLUDE_DIRS})
    endif()

    list(REMOVE_DUPLICATES G4_QT_INCLUDE_DIRS)
    list(REMOVE_ITEM G4_QT_INCLUDE_DIRS ${_cxx_compiler_dirs})

    set(G4_QT_CFLAGS )
    foreach(_dir ${G4_QT_INCLUDE_DIRS})
      set(G4_QT_CFLAGS "${G4_QT_CFLAGS} -I${_dir}")
    endforeach()

  else()
    set(G4_BUILTWITH_QT "no")
  endif()

  # - Wt
  if(GEANT4_USE_WT)
    set(G4_BUILTWITH_WT "yes")
    set(G4_WT_INCLUDE_DIRS ${Wt_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

    set(G4_WT_CFLAGS )
    foreach(_dir ${G4_WT_INCLUDE_DIRS})
      set(G4_WT_CFLAGS "${G4_WT_CFLAGS} -I${_dir}")
    endforeach()

  else()
    set(G4_BUILTWITH_WT "no")
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

  # If we have a module that uses X11, We have to play with the X11
  # paths to get a clean set suitable for inclusion
  if(G4_CONFIG_NEEDS_X11)
    set(_raw_x11_includes ${X11_INCLUDE_DIR})
    list(REMOVE_DUPLICATES _raw_x11_includes)
    list(REMOVE_ITEM _raw_x11_includes ${_cxx_compiler_dirs})
    set(G4_X11_CFLAGS )
    foreach(_p ${_raw_x11_includes})
      set(G4_X11_CFLAGS "-I${_p} ${G4_X11_CFLAGS}")
    endforeach()
  endif()

  # Configure the script
  # - BUILD TREE
  # Ouch, the include path will be LONG, but at least we always have
  # absolute paths...
  set(GEANT4_CONFIG_SELF_LOCATION "# BUILD TREE IS NON-RELOCATABLE")
  set(GEANT4_CONFIG_INSTALL_PREFIX "${PROJECT_BINARY_DIR}")
  set(GEANT4_CONFIG_INSTALL_EXECPREFIX \"\")
  # NB: this only works for *single* mode generators. With multimode
  # generators, which mode to use is not clear...
  set(GEANT4_CONFIG_LIBDIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

  get_property(__geant4_buildtree_include_dirs GLOBAL PROPERTY
    GEANT4_BUILDTREE_INCLUDE_DIRS)

  foreach(_dir ${__geant4_buildtree_include_dirs})
    set(GEANT4_CONFIG_INCLUDE_DIRS "${GEANT4_CONFIG_INCLUDE_DIRS} \\
    ${_dir}")
  endforeach()

  # - Data
  geant4_export_datasets(BUILD GEANT4_CONFIG_DATASET_DESCRIPTIONS)

  # Configure the build tree script
  # If we're on CMake 2.8 and above, we try to use file(COPY) to create an
  # executable script
  # Not sure if version check is o.k., but I'll be shocked if we ever see
  # a CMake 2.7 in the wild...
  if(${CMAKE_VERSION} VERSION_GREATER 2.7)
    configure_file(
      ${CMAKE_SOURCE_DIR}/cmake/Templates/geant4-config.in
      ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/geant4-config
      @ONLY
      )

    file(COPY
      ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/geant4-config
      DESTINATION ${PROJECT_BINARY_DIR}
      FILE_PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE
      )
  else()
    # Changing permissions is awkward, so just configure and document
    # that you have to do 'sh geant4-config' in this case.
    configure_file(
      ${CMAKE_SOURCE_DIR}/cmake/Templates/geant4-config.in
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

  # - Data
  geant4_export_datasets(INSTALL GEANT4_CONFIG_DATASET_DESCRIPTIONS)

  # Configure the install tree script
  configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/Templates/geant4-config.in
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

