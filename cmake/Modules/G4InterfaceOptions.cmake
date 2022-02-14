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

#-----------------------------------------------------------------------
# Select options first because these are interdependent
#
# - Coin3D/So{XT,Win,Qt}
option(GEANT4_USE_INVENTOR "Build Geant4 OpenInventor Xt/Win Visualization Driver" OFF)
geant4_add_feature(GEANT4_USE_INVENTOR "Build OpenInventor Xt/Win Driver")

option(GEANT4_USE_INVENTOR_QT "Build Geant4 OpenInventor Qt Visualization Driver" OFF)
geant4_add_feature(GEANT4_USE_INVENTOR_QT "Build OpenInventor Qt Driver")

# GEANT4_USE_INVENTOR is a dual-purpose variable - it marks use of Xt/Win driver for
# user, and marks enablement of G4OpenInventor internally. Check that user-set options
# are consistent before setting *internal* value
if(GEANT4_USE_INVENTOR AND GEANT4_USE_INVENTOR_QT)
  message(FATAL_ERROR "Only one of GEANT4_USE_INVENTOR and GEANT4_USE_INVENTOR_QT can be set ON, but both are ON.")
endif()

# User settings o.k., so mark that we're enabling G4OpenInventor
if(GEANT4_USE_INVENTOR_QT)
  set(GEANT4_USE_INVENTOR ON)
endif()

# - ToolsSG/Xt/X11/Qt/Win
# Though this option has platform differences, we configure here because it can
# trigger the Qt requirement
if(WIN32)
  set(_g4_toolssg_backends "OFF" "WIN32")
else()
  # TODO: Difference between X11, XT?
  set(_g4_toolssg_backends "OFF" "X11" "XT" "QT")
endif()

enum_option(GEANT4_USE_TOOLSSG
  DOC "Build Geant4 ToolsSG Visualization Driver with Backend"
  VALUES ${_g4_toolssg_backends}
)
geant4_add_feature(GEANT4_USE_TOOLSSG "Build ToolsSG Visualization Driver with backend '${GEANT4_USE_TOOLSSG}'")

foreach(_tsg_backend IN LISTS _g4_toolssg_backends)
  set(GEANT4_USE_TOOLSSG_${_tsg_backend}_GLES OFF)
  if(_tsg_backend STREQUAL GEANT4_USE_TOOLSSG)
    set(GEANT4_USE_TOOLSSG_${_tsg_backend}_GLES ON)
  endif()
endforeach()

# - Qt (may be required by Coin/ToolsSG driver)
option(GEANT4_USE_QT "Build Geant4 with Qt support" OFF)
if((GEANT4_USE_INVENTOR_QT OR GEANT4_USE_TOOLSSG_QT_GLES) AND NOT GEANT4_USE_QT)
  set(GEANT4_USE_QT ON CACHE BOOL "Build Geant4 with Qt support" FORCE)
  message(STATUS "Forcing GEANT4_USE_QT to ON, required by Inventor/ToolsSG Driver(s)")
endif()

geant4_add_feature(GEANT4_USE_QT "Build Geant4 with Qt support")

# - Vtk
option(GEANT4_USE_VTK "Build Geant4 with VTK visualisation" OFF)
mark_as_advanced(GEANT4_USE_VTK)
if(GEANT4_USE_VTK)
  find_package(VTK 8.2 REQUIRED)
endif()

geant4_add_feature(GEANT4_USE_VTK "Using VTK for visualisation (EXPERIMENTAL)")

# - Unix only
if(UNIX)
  # - Motif UI/Vis (may be required by Coin Vis driver)
  option(GEANT4_USE_XM "Build Geant4 with Motif (X11) support" OFF)
  if((GEANT4_USE_INVENTOR AND NOT GEANT4_USE_INVENTOR_QT)
      AND NOT GEANT4_USE_XM)
    set(GEANT4_USE_XM ON CACHE BOOL "Build Geant4 with Motif (X11) support" FORCE)
    message(STATUS "Forcing GEANT4_USE_XM to ON, required by Inventor driver")
  endif()

  geant4_add_feature(GEANT4_USE_XM "Build Geant4 with Xm Support")

  # OpenGL/X11 Vis Driver
  option(GEANT4_USE_OPENGL_X11 "Build Geant4 OpenGL driver with X11 support" OFF)
  geant4_add_feature(GEANT4_USE_OPENGL_X11 "Build Geant4 OpenGL driver with X11 support")

  # RayTracer driver with X11 support
  option(GEANT4_USE_RAYTRACER_X11 "Build RayTracer driver with X11 support" OFF)
  geant4_add_feature(GEANT4_USE_RAYTRACER_X11 "Build RayTracer driver with X11 support")
endif()

# Windows only
if(WIN32)
  option(GEANT4_USE_OPENGL_WIN32 "Build OpenGL driver with Win32 support")
  geant4_add_feature(GEANT4_USE_OPENGL_WIN32 "Build OpenGL driver with Win32 support")
endif()

#-----------------------------------------------------------------------
# Find dependencies
#
# Prefer Legacy GL when using CMake >= 3.10
set(OpenGL_GL_PREFERENCE "LEGACY")

# - Coin plus So{Xt,Qt,Win}
if(GEANT4_USE_INVENTOR)
  find_package(Coin 4.0.0 REQUIRED)
  geant4_save_package_variables(Inventor Coin_DIR)

  if(GEANT4_USE_INVENTOR_QT)
    find_package(SoQt 1.6.0 REQUIRED)
    geant4_save_package_variables(Inventor SoQt_DIR)
  else()
    if(UNIX)
      find_package(SoXt 1.4.0 REQUIRED)
      geant4_save_package_variables(Inventor SoXt_DIR)
      set(GEANT4_USE_INVENTOR_XT ON)
    elseif(WIN32)
      find_package(SoWin 1.4.0 REQUIRED)
      geant4_save_package_variables(Inventor SoWin_DIR)
      set(GEANT4_USE_INVENTOR_WIN ON)
    endif()
  endif()
endif()

# - Qt5
if(GEANT4_USE_QT)
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

  # Qt3D is only supported on 5.15 and above, but always on if available
  set(QT3D_MINIMUM_VERSION 5.15.0)
  set(GEANT4_USE_QT3D OFF)

  if(Qt5Core_VERSION GREATER_EQUAL QT3D_MINIMUM_VERSION)
    find_package(Qt53DCore ${QT3D_MINIMUM_VERSION} QUIET)
    find_package(Qt53DExtras ${QT3D_MINIMUM_VERSION} QUIET)
    find_package(Qt53DRender ${QT3D_MINIMUM_VERSION} QUIET)

    if(Qt53DCore_FOUND AND Qt53DExtras_FOUND AND Qt53DRender_FOUND)
      set(GEANT4_USE_QT3D ON)
      geant4_save_package_variables(Qt5 Qt53DCore_DIR Qt53DExtras_DIR Qt53DRender_DIR)
      geant4_add_feature(GEANT4_USE_QT3D "Build Geant4 Qt3D driver")
    else()
      set(_g4_qt53d_missing)
      if(NOT Qt53DCore_FOUND)
        list(APPEND _g4_qt53d_missing "Qt53DCore")
      endif()
      if(NOT Qt53DExtras_FOUND)
        list(APPEND _g4_qt53d_missing "Qt53DExtras")
      endif()
      if(NOT Qt53DRender_FOUND)
        list(APPEND _g4_qt53d_missing "Qt53DRender")
      endif()

      message(STATUS "Disabling Geant4 Qt3D driver, missing Qt5 packages: ${_g4_qt53d_missing}")
    endif()
  endif()
  # Variables for export to GNUmake
  execute_process(COMMAND ${QT_QMAKE_EXECUTABLE} -query QT_INSTALL_PREFIX OUTPUT_VARIABLE G4QTHOME OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${QT_QMAKE_EXECUTABLE} -query QT_INSTALL_LIBS OUTPUT_VARIABLE G4QTLIBPATH OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

# - OpenGL
if(GEANT4_USE_INVENTOR
   OR GEANT4_USE_QT
   OR GEANT4_USE_OPENGL_X11
   OR GEANT4_USE_XM
   OR GEANT4_USE_TOOLSSG
   OR GEANT4_USE_OPENGL_WIN32)
  find_package(OpenGL REQUIRED)
  geant4_save_package_variables(OpenGL OPENGL_INCLUDE_DIR OPENGL_gl_LIBRARY)

  # X11 drivers on macOS need XQuartzGL
  if(APPLE AND
     (GEANT4_USE_INVENTOR_XT
      OR GEANT4_USE_OPENGL_X11
      OR GEANT4_USE_XM
      OR GEANT4_USE_TOOLSSG_X11_GLES
      OR GEANT4_USE_TOOLSSG_XT_GLES))
    find_package(XQuartzGL REQUIRED)
    geant4_save_package_variables(XQuartzGL XQuartzGL_INCLUDE_DIR XQuartzGL_gl_LIBRARY)
  endif()

  # Enable G4OpenGL driver
  if(GEANT4_USE_INVENTOR
     OR GEANT4_USE_QT
     OR GEANT4_USE_OPENGL_X11
     OR GEANT4_USE_OPENGL_WIN32)
    set(GEANT4_USE_OPENGL ON)
  endif()
endif()

if(UNIX)
  # - X11
  if(GEANT4_USE_OPENGL_X11
     OR GEANT4_USE_RAYTRACER_X11
     OR GEANT4_USE_XM
     OR GEANT4_USE_INVENTOR_XT
     OR GEANT4_USE_TOOLSSG_X11_GLES
     OR GEANT4_USE_TOOLSSG_XT_GLES)
    find_package(X11 REQUIRED)
    include("${CMAKE_CURRENT_LIST_DIR}/G4X11Shim.cmake")

    # Check for additional required X11 libraries
    if(NOT X11_Xmu_FOUND)
      message(FATAL_ERROR "could not find X11 Xmu library and/or headers")
    endif()

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
      X11_Xt_LIB)
  endif()

  # - Motif
  if(GEANT4_USE_XM)
    find_package(Motif REQUIRED)
    include("${CMAKE_CURRENT_LIST_DIR}/G4MotifShim.cmake")
    geant4_save_package_variables(Motif MOTIF_INCLUDE_DIR MOTIF_LIBRARIES)
  endif()
endif()

