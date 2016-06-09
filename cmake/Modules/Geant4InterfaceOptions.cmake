# - Configure options for Geant4 UI/Vis requiring Third Party Libraries
# All core Interfaces not requiring third party support are built 
# automatically with platform differences taken into account. 
# These are:
#
# UI Category : Always built on supported platforms
#  G4UIterminal : UNIX, WIN32
#    + - G4UItcsh : UNIX only.
#    + - G4UIcsh  : UNIX, WIN32
#  G4UIWin32    : WIN32
#
# UI category : Built on supported platforms when Third Party 
# library(libraries) available
#  G4UIQt  : UNIX, WIN32 (Requires : Qt4)
#  G4UIXm  : UNIX (Requires : Qt4)
#  G4UIXaw : DEPRECATED
#
# Vis Category: Always built on supported platforms
#  DAWNFILE : UNIX, Win32
#  ...
#  gMocren
#  HepRep
#  RayTracer
#  Tree
#  VRML
#
# Vis Category : Built on supported platforms when Third Party 
# library(libraries) available
#  OpenGL : UNIX, WIN32
#  OpenInventor : UNIX WIN32
#
# Rather than separate the UI and visualization options, we simply ask if
# the user wants to build support for a particular tool, e.g. Qt.
#
# Specific UI/Vis options are also handled here.
#


#-----------------------------------------------------------------------
# Configure OpenInventor support
#
option(GEANT4_USE_INVENTOR "Build Geant4 OpenInventor Visualization Driver" OFF)

if(GEANT4_USE_INVENTOR)
  # Always need Inventor and OpenGL
  find_package(Inventor REQUIRED)
  find_package(OpenGL REQUIRED)

  # On UNIX, we also require Xm and X11 (inc. Xpm)...
  if(UNIX)
    if(NOT INVENTOR_SOXT_LIBRARY OR NOT INVENTOR_SOXT_INCLUDE_DIR)
      message(FATAL_ERROR "Could not find SoXt library and/or headers")
    endif()
    set(INVENTOR_INCLUDE_DIR "${INVENTOR_INCLUDE_DIR}" "${INVENTOR_SOXT_INCLUDE_DIR}")
    find_package(Motif REQUIRED)
    find_package(X11 REQUIRED)

    if(NOT X11_Xpm_FOUND)
      message(FATAL_ERROR "Could not find X11 Xpm headers and/or library (Required by GEANT4_USE_INVENTOR)")
    endif()
  endif()

  GEANT4_ADD_FEATURE(
    GEANT4_USE_INVENTOR "Build OpenInventor Driver"
    )
endif()


#-----------------------------------------------------------------------
# Configure Qt Support if needed (CROSSPLATFORM).
#
option(GEANT4_USE_QT "Build Geant4 with Qt4 support" OFF)

if(GEANT4_USE_QT)
  # Find and configure Qt and OpenGL - require 4
  # This is fine on Mac OS X because Qt will use Framework GL.
  # On WIN32 only, set QT_USE_IMPORTED_TARGETS. The Qt4 module will 
  # otherwise set QT_LIBRARIES using the 'optimized A; debug Ad' 
  # technique. CMake will then complain because when building DLLs we reset 
  # LINK_INTERFACE_LIBRARIES and this cannot be passed a link of this form.
  # This means we have to recreate the imported targets later though...
  if(WIN32)
    set(QT_USE_IMPORTED_TARGETS ON)
  endif()

  find_package(Qt4 REQUIRED COMPONENTS QtCore QtGui QtOpenGL)
  find_package(OpenGL REQUIRED)

  # OpenGL part of Qt is in OpenGL component so mark the need to
  # add OpenGL.
  set(GEANT4_USE_OPENGL ON)
  GEANT4_ADD_FEATURE(GEANT4_USE_QT "Build Geant4 with Qt support")
endif()


#-----------------------------------------------------------------------
# UNIX PLATFORMS ONLY
#-----------------------------------------------------------------------
if(UNIX)
  # - Support for Client/Server DAWN driver
  # mark as advanced because user should know what they're doing to use this
  option(GEANT4_USE_NETWORKDAWN "Build Dawn driver with Client/Server support" OFF)
  #
  # Possible headers checks for needed network parts?
  #
  GEANT4_ADD_FEATURE(GEANT4_USE_NETWORKDAWN "Build Dawn driver with Client/Server support")
  mark_as_advanced(GEANT4_USE_NETWORKDAWN)

  # - Support for Client/Server VRML driver
  # mark as advanced because user should know what they're doing to use this
  option(GEANT4_USE_NETWORKVRML "Build VRML driver with Client/Server support" OFF)
  #
  # Possible header checks for needed network parts?
  #
  GEANT4_ADD_FEATURE(GEANT4_USE_NETWORKVRML "Build VRML driver with Client/Server support")
  mark_as_advanced(GEANT4_USE_NETWORKVRML)

  # - Configure Xm support if requested
  option(GEANT4_USE_XM "Build Geant4 with Motif (X11) support" OFF)

  # - Configure OpenGL X11 support if requested
  option(GEANT4_USE_OPENGL_X11 "Build Geant4 OpenGL driver with X11 support"
    OFF)

  # - X11 Support for RayTracer driver
  option(GEANT4_USE_RAYTRACER_X11 "Build RayTracer driver with X11 support" OFF)

  # -- Now configure needed X11, Motif and OpenGL packages
  # We also have to be concerned with Inventor on UNIX, because
  # that also needs the Xt configuration...
  # - X11
  if(GEANT4_USE_XM OR GEANT4_USE_OPENGL_X11 OR GEANT4_USE_RAYTRACER_X11 OR GEANT4_USE_INVENTOR)
    find_package(X11 REQUIRED)

    # We also require Xmu, which isn't found by default
    # Just the search in here, copying pattern from 
    # FindX11. We don't add it to X11 libraries because it's only
    # needed in the OpenGL driver.
    #
    set(CMAKE_FIND_FRAMEWORK_SAVE ${CMAKE_FIND_FRAMEWORK})
    set(CMAKE_FIND_FRAMEWORK NEVER)

    set(X11_INC_SEARCH_PATH
        /usr/pkg/xorg/include
        /usr/X11R6/include 
        /usr/X11R7/include 
        /usr/include/X11
        /usr/openwin/include 
        /usr/openwin/share/include 
        /opt/graphics/OpenGL/include
        )

    set(X11_LIB_SEARCH_PATH
        /usr/pkg/xorg/lib
        /usr/X11R6/lib
        /usr/X11R7/lib
        /usr/openwin/lib 
        )

    find_path(X11_Xmu_INCLUDE_PATH X11/Xmu/Xmu.h ${X11_INC_SEARCH_PATH})
    find_library(X11_Xmu_LIBRARY Xmu ${X11_SEARCH_PATH})
    if(X11_Xmu_LIBRARY-NOTFOUND OR X11_Xmu_INCLUDE_PATH-NOTFOUND)
      message(FATAL_ERROR "could not find X11 Xmu library and/or headers")
    endif()

    mark_as_advanced(X11_Xmu_INCLUDE_PATH X11_Xmu_LIBRARY)
    set(CMAKE_FIND_FRAMEWORK ${CMAKE_FIND_FRAMEWORK_SAVE})
    # END of finding Xmu

    # Check for additional required X11 libraries
    # - Xt
    if(NOT X11_Xt_FOUND)
      message(FATAL_ERROR "could not find X11 Xt library and/or headers")
    endif()

    # Create final list of X11 libraries needed for Geant4
    # - Motif
    if(APPLE)
      if(GEANT4_USE_XM OR GEANT4_USE_INVENTOR)
        set(X11_LIBRARIES ${X11_LIBRARIES} ${X11_Xt_LIB})
      endif()
    endif()

    # - If we got to this point, RayTracer is o.k., so add the feature
    GEANT4_ADD_FEATURE(GEANT4_USE_RAYTRACER_X11 "Build RayTracer driver with X11 support")
  endif()

  #- OpenGL : we also do this for Inventor so that we pick up Mac X11 GL...
  if(GEANT4_USE_INVENTOR OR GEANT4_USE_XM OR GEANT4_USE_OPENGL_X11)
    # - Find native GL first
    find_package(OpenGL REQUIRED)

    # - If we're on Apple, also find the X11 GL libraries and append them.
    if(APPLE)
      # - This is for X11 GL drivers, so we DON'T want Framework!
      set(CMAKE_FIND_FRAMEWORK_SAVE ${CMAKE_FIND_FRAMEWORK})
      set(CMAKE_FIND_FRAMEWORK NEVER)

      find_path(OPENGL_X11_INCLUDE_DIR GL/gl.h
        PATHS /usr/X11R6/include
        NO_DEFAULT_PATH
        )

      find_library(OPENGL_X11_gl_LIBRARY GL
        PATHS /usr/X11R6/lib
        NO_DEFAULT_PATH
        )

      find_library(OPENGL_X11_glu_LIBRARY GLU
        PATHS /usr/X11R6/lib
        NO_DEFAULT_PATH
        )

      mark_as_advanced(
        OPENGL_X11_INCLUDE_DIR 
        OPENGL_X11_gl_LIBRARY 
        OPENGL_X11_glu_LIBRARY
        )

      if(NOT OPENGL_X11_INCLUDE_DIR)
        message(FATAL_ERROR "Apple X11 OpenGL GL/gl.h NOTFOUND")
      endif()

      if(NOT OPENGL_X11_gl_LIBRARY)
        message(FATAL_ERROR "Apple X11 OpenGL GL Library NOTFOUND")
      endif()

      if(NOT OPENGL_X11_glu_LIBRARY)
        message(FATAL_ERROR "Apple X11 OpenGL GLU Library NOTFOUND")
      endif()

      # Append the X11 GL libraries to the native paths
      set(OPENGL_INCLUDE_DIR ${OPENGL_INCLUDE_DIR} ${OPENGL_X11_INCLUDE_DIR})
      set(OPENGL_LIBRARIES ${OPENGL_LIBRARIES} 
        ${OPENGL_X11_gl_LIBRARY}
        ${OPENGL_X11_glu_LIBRARY}
        )

      # Restore framework search
      set(CMAKE_FIND_FRAMEWORK ${CMAKE_FIND_FRAMEWORK_SAVE})
    endif()

    # Add feature for X11 GL
    GEANT4_ADD_FEATURE(GEANT4_USE_OPENGL_X11 "Build Geant4 OpenGL driver with X11 support")
  endif()

  # - Add need for OpenGL in X11 case
  if(GEANT4_USE_OPENGL_X11 OR GEANT4_USE_XM)
    set(GEANT4_USE_OPENGL ON)
  endif()

  # - Motif last of all, then can add the feature
  if(GEANT4_USE_XM)
    find_package(Motif REQUIRED)
    GEANT4_ADD_FEATURE(GEANT4_USE_XM "Build Geant4 with Xm Support")
  endif()  
endif()

#-----------------------------------------------------------------------
# WINDOWS PLATFORMS ONLY
#-----------------------------------------------------------------------
if(WIN32)
  # - OpenGL Win32...
  option(GEANT4_USE_OPENGL_WIN32 "Build OpenGL driver with Win32 support")

  if(GEANT4_USE_OPENGL_WIN32)
    # Just need OpenGL and on Windows, this should be easy...
    find_package(OpenGL REQUIRED)

    # This is part of the G4OpenGL component, so mark it as needed
    set(GEANT4_USE_OPENGL ON)
    GEANT4_ADD_FEATURE(
      GEANT4_USE_OPENGL_WIN32 "Build OpenGL driver with Win32 support"
      )
  endif()
endif()

# - and we're done...

