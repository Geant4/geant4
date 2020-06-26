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
#  G4UIQt  : UNIX, WIN32 (Requires : Qt5)
#  G4UIWt  : Web (Requires : Wt) DEPRECATED
#  G4UIXm  : UNIX
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

# Prefer Legacy GL when using CMake >= 3.10
set(OpenGL_GL_PREFERENCE "LEGACY")

#-----------------------------------------------------------------------
# Configure OpenInventor support
#
option(GEANT4_USE_INVENTOR "Build Geant4 OpenInventor Visualization Driver" OFF)
option(GEANT4_USE_INVENTOR_QT "Build Geant4 OpenInventor Qt Visualization Driver" OFF)

if(GEANT4_USE_INVENTOR_QT)
  set(GEANT4_USE_INVENTOR ON)
  set(GEANT4_USE_QT ON)
endif()

if(GEANT4_USE_INVENTOR)
  # Always need Inventor and OpenGL
  find_package(Inventor REQUIRED)
  find_package(OpenGL REQUIRED)

  # On UNIX, we also require Xm and X11 (inc. Xpm)...
  if(UNIX)
    if (GEANT4_USE_INVENTOR_QT)
      if (NOT INVENTOR_SOQT_LIBRARY OR NOT INVENTOR_SOQT_INCLUDE_DIR)
        message(FATAL_ERROR "Could not find SoQt library and/or headers")
      endif()
      set(INVENTOR_INCLUDE_DIR "${INVENTOR_INCLUDE_DIR}" "${INVENTOR_SOQT_INCLUDE_DIR}")
    else()
      if(NOT INVENTOR_SOXT_LIBRARY OR NOT INVENTOR_SOXT_INCLUDE_DIR)
        message(FATAL_ERROR "Could not find SoXt library and/or headers")
      endif()
      set(INVENTOR_INCLUDE_DIR "${INVENTOR_INCLUDE_DIR}" "${INVENTOR_SOXT_INCLUDE_DIR}")
      set(GEANT4_USE_XM ON)
      find_package(Motif REQUIRED)
    endif()

    find_package(X11 REQUIRED)

    if(NOT X11_Xpm_FOUND)
      message(FATAL_ERROR "Could not find X11 Xpm headers and/or library (Required by GEANT4_USE_INVENTOR)")
    endif()
  endif()

  geant4_add_feature(GEANT4_USE_INVENTOR "Build OpenInventor Driver")
endif()


#-----------------------------------------------------------------------
# Configure Qt Support if needed (CROSSPLATFORM).
#
option(GEANT4_USE_QT "Build Geant4 with Qt support" OFF)

if(GEANT4_USE_QT)
  # Find and configure Qt and OpenGL
  find_package(Qt5Core REQUIRED)
  find_package(Qt5Gui REQUIRED)
  find_package(Qt5Widgets REQUIRED)
  find_package(Qt5OpenGL REQUIRED)
  find_package(Qt5PrintSupport REQUIRED)

  set(QT_QTCORE_LIBRARY Qt5::Core Qt5::PrintSupport)
  set(QT_QTGUI_LIBRARY Qt5::Gui Qt5::Widgets)
  set(QT_OPENGL_LIBRARY Qt5::Gui Qt5::OpenGL)
  set(QT_LIBRARIES Qt5::OpenGL Qt5::Gui Qt5::PrintSupport Qt5::Widgets)
  get_target_property(QT_QMAKE_EXECUTABLE ${Qt5Core_QMAKE_EXECUTABLE} IMPORTED_LOCATION)
  geant4_save_package_variables(Qt5 Qt5Core_DIR Qt5Gui_DIR Qt5Widgets_DIR Qt5OpenGL_DIR Qt5PrintSupport_DIR)

  find_package(OpenGL REQUIRED)
  geant4_save_package_variables(OpenGL OPENGL_INCLUDE_DIR OPENGL_gl_LIBRARY)

  # Variables for export to GNUmake
  execute_process(COMMAND ${QT_QMAKE_EXECUTABLE} -query QT_INSTALL_PREFIX OUTPUT_VARIABLE G4QTHOME OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${QT_QMAKE_EXECUTABLE} -query QT_INSTALL_LIBS OUTPUT_VARIABLE G4QTLIBPATH OUTPUT_STRIP_TRAILING_WHITESPACE)

  # OpenGL part of Qt is in OpenGL component so mark the need to
  # add OpenGL.
  set(GEANT4_USE_OPENGL ON)
  geant4_add_feature(GEANT4_USE_QT "Build Geant4 with Qt support")
endif()

#-----------------------------------------------------------------------
# Configure Wt Support
#
if(GEANT4_USE_WT)
  message(WARNING "Support for Wt in Geant4 is no longer available and may be removed in future releases")
endif()
set(GEANT4_USE_WT OFF)

#option(GEANT4_USE_WT "Build Geant4 with Wt4 support" OFF)
#mark_as_advanced(GEANT4_USE_WT)
#
#if(GEANT4_USE_WT)
#  # Find and configure Wt and OpenGL
#  find_package(Wt REQUIRED)
#  find_package(OpenGL REQUIRED)
#
#  geant4_save_package_variables(Wt Wt_INCLUDE_DIR Wt_LIBRARY_RELEASE Wt_LIBRARY_DEBUG Wt_HTTP_LIBRARY_RELEASE Wt_HTTP_LIBRARY_DEBUG)
#  geant4_save_package_variables(OpenGL OPENGL_INCLUDE_DIR OPENGL_gl_LIBRARY)
#
#  set(WT_DEFINITIONS "-DQT_NO_KEYWORDS")
#
#  # WebGL part of Wt is in OpenGL component so mark the need to
#  # add OpenGL.
#  set(GEANT4_USE_OPENGL ON)
#  geant4_add_feature(GEANT4_USE_WT "Build Geant4 with Wt support")
#endif()


#-----------------------------------------------------------------------
# UNIX PLATFORMS ONLY
#-----------------------------------------------------------------------
if(UNIX)
  # - Support for Client/Server DAWN driver
  # mark as advanced because user should know what they're doing here
  option(GEANT4_USE_NETWORKDAWN "Build Dawn driver with Client/Server support" OFF)
  #
  # Possible headers checks for needed network parts?
  #
  geant4_add_feature(GEANT4_USE_NETWORKDAWN "Build Dawn driver with Client/Server support")
  mark_as_advanced(GEANT4_USE_NETWORKDAWN)

  # - Support for Client/Server VRML driver
  # mark as advanced because user should know what they're doing to use this
  option(GEANT4_USE_NETWORKVRML "Build VRML driver with Client/Server support" OFF)
  geant4_add_feature(GEANT4_USE_NETWORKVRML "Build VRML driver with Client/Server support")
  mark_as_advanced(GEANT4_USE_NETWORKVRML)

  # - Configure Xm support if requested
  option(GEANT4_USE_XM "Build Geant4 with Motif (X11) support" OFF)

  # - Configure OpenGL X11 support if requested
  option(GEANT4_USE_OPENGL_X11 "Build Geant4 OpenGL driver with X11 support" OFF)

  # - X11 Support for RayTracer driver
  option(GEANT4_USE_RAYTRACER_X11 "Build RayTracer driver with X11 support" OFF)

  # -- Now configure needed X11, Motif and OpenGL packages
  # We also have to be concerned with Inventor on UNIX, because
  # that also needs the Xt configuration...
  # - X11
  if(GEANT4_USE_XM OR GEANT4_USE_OPENGL_X11 OR GEANT4_USE_RAYTRACER_X11 OR GEANT4_USE_INVENTOR)
    find_package(X11 REQUIRED)
    include("${CMAKE_CURRENT_LIST_DIR}/G4X11Shim.cmake")

    # Check for additional required X11 libraries
    # - Xmu
    if(NOT X11_Xmu_FOUND)
      message(FATAL_ERROR "could not find X11 Xmu library and/or headers")
    endif()

    # - Xt
    if(NOT X11_Xt_FOUND)
      message(FATAL_ERROR "could not find X11 Xt library and/or headers")
    endif()

    geant4_save_package_variables(X11
      X11_X11_INCLUDE_PATH
      X11_X11_LIB
      X11_ICE_INCLUDE_PATH
      X11_ICE_LIB
      X11_SM_INCLUDE_PATH
      X11_SM_LIB
      X11_Xext_INCLUDE_PATH
      X11_Xext_LIB
      X11_Xmu_INCLUDE_PATH
      X11_Xmu_LIB
      X11_Xt_INCLUDE_PATH
      X11_Xt_LIB
    )

    # - If we got to this point, RayTracer is o.k., so add the feature
    geant4_add_feature(GEANT4_USE_RAYTRACER_X11 "Build RayTracer driver with X11 support")
  endif()

  #- OpenGL : we also do this for Inventor so that we pick up Mac X11 GL...
  if(GEANT4_USE_INVENTOR OR GEANT4_USE_XM OR GEANT4_USE_OPENGL_X11)
    # - Find native GL first
    find_package(OpenGL REQUIRED)
    geant4_save_package_variables(OpenGL OPENGL_INCLUDE_DIR OPENGL_gl_LIBRARY)

    # - If we're on Apple, also find the XQuartz GL libraries and append t.
    if(APPLE)
      find_package(XQuartzGL REQUIRED)
      geant4_save_package_variables(XQuartzGL XQuartzGL_INCLUDE_DIR XQuartzGL_gl_LIBRARY)
    endif()

    # Add feature for X11 GL
    geant4_add_feature(GEANT4_USE_OPENGL_X11 "Build Geant4 OpenGL driver with X11 support")
  endif()

  # - Add need for OpenGL in X11 case
  if(GEANT4_USE_OPENGL_X11 OR GEANT4_USE_XM)
    set(GEANT4_USE_OPENGL ON)
  endif()

  # - Motif last of all, then can add the feature
  if(GEANT4_USE_XM)
    find_package(Motif REQUIRED)
    include("${CMAKE_CURRENT_LIST_DIR}/G4MotifShim.cmake")
    geant4_save_package_variables(Motif MOTIF_INCLUDE_DIR MOTIF_LIBRARIES)
    geant4_add_feature(GEANT4_USE_XM "Build Geant4 with Xm Support")
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
    geant4_add_feature(GEANT4_USE_OPENGL_WIN32 "Build OpenGL driver with Win32 support")
    geant4_save_package_variables(OpenGL OPENGL_INCLUDE_DIR OPENGL_gl_LIBRARY)
  endif()
endif()

# - and we're done...

