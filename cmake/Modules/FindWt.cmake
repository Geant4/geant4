# Find Wt includes and libraries
#
# This script sets the following variables:
#
#  Wt_INCLUDE_DIR
#  Wt_LIBRARIES  - Release libraries
#  Wt_FOUND  - True if release libraries found
#  WT_FOUND  - True if release libraries found
#  Wt_DEBUG_LIBRARIES  - Debug libraries
#  Wt_DEBUG_FOUND  - True if debug libraries found
#
# To direct the script to a particular Wt installation, use the
# standard cmake variables CMAKE_INCLUDE_PATH and CMAKE_LIBRARY_PATH
#
# To use this script to find Wt, when using the new style for include files:
#   #include <Wt/WLineEdit>
#   #include <Wt/Ext/LineEdit>
#   #include <Wt/Chart/WPieChart>
#
# include the following CMake snippet in your project:
#
#  FIND_PACKAGE( Wt REQUIRED )
#  INCLUDE_DIRECTORIES( ${Wt_INCLUDE_DIR} )
#  TARGET_LINK_LIBRARIES( yourexe
#    ${Wt_DEBUG_LIBRARY}        # or {Wt_LIBRARY}
#    ${Wt_HTTP_DEBUG_LIBRARY}   # or {Wt_HTTP_LIBRARY}
#    ${Wt_EXT_DEBUG_LIBRARY}    # or {Wt_EXT_LIBRARY}
#  )
#
# To use this script to find Wt, when using the old include style:
#   #include <WLineEdit>
#   #include <Ext/LineEdit>
#   #include <Chart/WPieChart>
# style of include files, change the INCLUDE_DIRECTORIES statement to:
#   INCLUDE_DIRECTORIES( ${Wt_INCLUDE_DIR} ${Wt_INCLUDE_DIR}/Wt )
#
#
#
#
# Copyright (c) 2007, Pau Garcia i Quiles, <pgquiles@elpauer.org>
#
# Modifications Copyright (c) 2013, Ben Morgan, <Ben.Morgan@warwick.ac.uk>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.


find_path(Wt_INCLUDE_DIR NAMES Wt/WObject PATHS ENV PATH PATH_SUFFIXES include wt )

include(SelectLibraryConfigurations)

find_library( Wt_LIBRARY_RELEASE NAMES wt PATHS PATH PATH_SUFFIXES lib lib-release lib_release )
find_library( Wt_LIBRARY_DEBUG NAMES wtd PATHS PATH PATH_SUFFIXES lib libd lib-debug lib_debug HINTS /usr/lib/debug/usr/lib)
select_library_configurations(Wt)

find_library( Wt_EXT_LIBRARY_RELEASE NAMES wtext PATHS PATH PATH_SUFFIXES lib lib-release lib_release )
find_library( Wt_EXT_LIBRARY_DEBUG NAMES wtextd PATHS PATH PATH_SUFFIXES lib libd lib-debug lib_debug HINTS /usr/lib/debug/usr/lib)
select_library_configurations(Wt_EXT)

find_library( Wt_HTTP_LIBRARY_RELEASE NAMES wthttp PATHS PATH PATH_SUFFIXES lib lib-release lib_release )
find_library( Wt_HTTP_LIBRARY_DEBUG NAMES wthttpd PATHS PATH PATH_SUFFIXES lib libd lib-debug lib_debug HINTS /usr/lib/debug/usr/lib)
select_library_configurations(Wt_HTTP)

find_library( Wt_FCGI_LIBRARY_RELEASE NAMES wtfcgi PATHS PATH PATH_SUFFIXES lib lib-release lib_release )
find_library( Wt_FCGI_LIBRARY_DEBUG NAMES wtfcgid PATHS PATH PATH_SUFFIXES lib libd lib-debug lib_debug HINTS /usr/lib/debug/usr/lib)
select_library_configurations(Wt_FCGI)

find_library( Wt_DBO_LIBRARY_RELEASE NAMES wtdbo PATHS PATH PATH_SUFFIXES lib lib-release lib_release )
find_library( Wt_DBO_LIBRARY_DEBUG NAMES wtdbod PATHS PATH PATH_SUFFIXES lib lib-debug lib_debug HINTS /usr/lib/debug/usr/lib)
select_library_configurations(Wt_DBO)

find_library( Wt_DBOSQLITE3_LIBRARY_RELEASE NAMES wtdbosqlite3 PATHS PATH PATH_SUFFIXES lib lib-release lib_release )
find_library( Wt_DBOSQLITE3_LIBRARY_DEBUG NAMES wtdbosqlite3d PATHS PATH PATH_SUFFIXES lib lib-debug lib_debug HINTS /usr/lib/debug/usr/lib)
select_library_configurations(Wt_DBOSQLITE3)

find_library( Wt_DBOPOSTGRES_LIBRARY_RELEASE NAMES wtdbopostgres PATHS PATH PATH_SUFFIXES lib lib-release lib_release )
find_library( Wt_DBOPOSTGRES_LIBRARY_DEBUG NAMES wtdbopostgresd PATHS PATH PATH_SUFFIXES lib lib-debug lib_debug HINTS /usr/lib/debug/usr/lib)
select_library_configurations(Wt_DBOPOSTGRES)
        
# Need at least Boost signals...
find_package(Boost REQUIRED signals)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Wt REQUIRED_VARS Wt_INCLUDE_DIR Wt_LIBRARY)

# - Create imported targets as needed by Geant4
if(Wt_FOUND)
  # Main library
  if(NOT TARGET Wt::Wt)
    add_library(Wt::Wt UNKNOWN IMPORTED)
    set_target_properties(Wt::Wt PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${Wt_INCLUDE_DIR}"
      INTERFACE_LINK_LIBRARIES "Boost::signals;Boost::headers")

    if(Wt_LIBRARY_RELEASE)
      set_property(TARGET Wt::Wt APPEND PROPERTY
        IMPORTED_CONFIGURATIONS RELEASE)
      set_target_properties(Wt::Wt PROPERTIES
        IMPORTED_LOCATION_RELEASE "${Wt_LIBRARY_RELEASE}")
    endif()

    if(Wt_LIBRARY_DEBUG)
      set_property(TARGET Wt::Wt APPEND PROPERTY
        IMPORTED_CONFIGURATIONS DEBUG)
      set_target_properties(Wt::Wt PROPERTIES
        IMPORTED_LOCATION_DEBUG "${Wt_LIBRARY_DEBUG}")
    endif()

    if(NOT Wt_LIBRARY_RELEASE AND NOT Wt_LIBRARY_DEBUG)
      set_property(TARGET Wt::Wt APPEND PROPERTY
        IMPORTED_LOCATION "${Wt_LIBRARY}")
     endif()
  endif()


  # Components
  foreach(__wt_target HTTP)
    if(NOT TARGET Wt::${__wt_target})
      add_library(Wt::${__wt_target} UNKNOWN IMPORTED)
      set_target_properties(Wt::${__wt_target} PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${Wt_INCLUDE_DIR}")

      if(Wt_${__wt_target}_LIBRARY_RELEASE)
        set_property(TARGET Wt::${__wt_target} APPEND PROPERTY
          IMPORTED_CONFIGURATIONS RELEASE)
        set_target_properties(Wt::${__wt_target} PROPERTIES
          IMPORTED_LOCATION_RELEASE "${Wt_${__wt_target}_LIBRARY_RELEASE}")
      endif()

      if(Wt_${__wt_target}_LIBRARY_DEBUG)
        set_property(TARGET Wt::${__wt_target} APPEND PROPERTY
          IMPORTED_CONFIGURATIONS DEBUG)
        set_target_properties(Wt::${__wt_target} PROPERTIES
          IMPORTED_LOCATION_DEBUG "${Wt_${__wt_target}_LIBRARY_DEBUG}")
      endif()

      if(NOT Wt_${__wt_target}_LIBRARY_RELEASE AND NOT Wt_${__wt_target}_LIBRARY_DEBUG)
        set_property(TARGET Wt::${__wt_target} APPEND PROPERTY
          IMPORTED_LOCATION "${Wt_${__wt_target}_LIBRARY}")
      endif()
    endif()
  endforeach()
endif()

