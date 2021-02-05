#[=======================================================================[.rst:
FindXQuartzGL
-------------

FindModule for XQuartz/Homebrew/MacPorts implementation of OpenGL/GLU. Specific
to Geant4 to allow use of XQuartz/Homebrew/MacPorts only.

Use of the module on non-macOS systems will result in a fatal error

IMPORTED Targets
^^^^^^^^^^^^^^^^

This module defines the :prop_tgt:`IMPORTED` targets:

``XQuartz::GL``
 Defined to the XQuartz GL library
``XQuartz::GLU``
 Define to the XQuartz GLU library

Result Variables
^^^^^^^^^^^^^^^^

This module sets the following variables:

``XQuartzGL_FOUND``
 True, if the XQuartz GL libraries were located

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set

``XQuartzGL_INCLUDE_DIR``
 Path to the XQuartz GL include directory
``XQuartzGL_gl_LIBRARY``
 Path to the XQuartz GL library
``XQuartzGL_glu_LIBRARY``
 Path to the XQuartz GLU library

#]=======================================================================]

# Just don't run if we're on macOS
if(NOT APPLE)
  message(FATAL_ERROR "FindXQuartzGL is only for use on macOS platforms")
endif()

# - This is for X11 GL drivers, so we DON'T want Framework!
set(CMAKE_FIND_FRAMEWORK_SAVE ${CMAKE_FIND_FRAMEWORK})
set(CMAKE_FIND_FRAMEWORK NEVER)

find_path(XQuartzGL_INCLUDE_DIR GL/gl.h
  PATHS /usr/X11R6/include /opt/X11/include /usr/local/include /opt/local/include
  NO_DEFAULT_PATH
  )

find_library(XQuartzGL_gl_LIBRARY GL
  PATHS /usr/X11R6/lib /opt/X11/lib /usr/local/lib /opt/local/lib
  NO_DEFAULT_PATH
  )

find_library(XQuartzGL_glu_LIBRARY GLU
  PATHS /usr/X11R6/lib /opt/X11/lib /usr/local/lib /opt/local/lib
  NO_DEFAULT_PATH
  )

set(CMAKE_FIND_FRAMEWORK ${CMAKE_FIND_FRAMEWORK_SAVE})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(XQuartzGL
  FOUND_VAR
    XQuartzGL_FOUND
  REQUIRED_VARS
    XQuartzGL_INCLUDE_DIR
    XQuartzGL_gl_LIBRARY
  )

mark_as_advanced(XQuartzGL_INCLUDE_DIR XQuartzGL_gl_LIBRARY XQuartzGL_glu_LIBRARY)

if(XQuartzGL_FOUND)
  if(NOT TARGET XQuartzGL::GL)
    add_library(XQuartzGL::GL UNKNOWN IMPORTED)
    set_target_properties(XQuartzGL::GL PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${XQuartzGL_INCLUDE_DIR}"
      IMPORTED_LOCATION "${XQuartzGL_gl_LIBRARY}"
    )
  endif()

  if(NOT TARGET XQuartzGL::GLU AND XQuartzGL_glu_LIBRARY)
    add_library(XQuartzGL::GLU UNKNOWN IMPORTED)
    set_target_properties(XQuartzGL::GLU PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${XQuartzGL_INCLUDE_DIR}"
      INTERFACE_LINK_LIBRARIES XQuartzGL::GL
      IMPORTED_LOCATION "${XQuartzGL_glu_LIBRARY}"
    )
  endif()
endif()
