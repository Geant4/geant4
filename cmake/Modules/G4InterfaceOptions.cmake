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

# - Qt
# May be required by Coin driver
# Selection also enables ToolsSG driver
option(GEANT4_USE_QT "Build Geant4 with Qt support" OFF)
if(GEANT4_USE_INVENTOR_QT AND NOT GEANT4_USE_QT)
  set(GEANT4_USE_QT ON CACHE BOOL "Build Geant4 with Qt support" FORCE)
  message(STATUS "Forcing GEANT4_USE_QT to ON, required by selection of GEANT4_USE_INVENTOR_QT as ON")
endif()
set(GEANT4_USE_TOOLSSG_QT_GLES ${GEANT4_USE_QT})

# - Vtk
option(GEANT4_USE_VTK "Build Geant4 with VTK visualisation" OFF)
mark_as_advanced(GEANT4_USE_VTK)
if(GEANT4_USE_VTK)
  find_package(VTK 8.2 REQUIRED)
endif()
geant4_add_feature(GEANT4_USE_VTK "Using VTK for visualisation (EXPERIMENTAL)")

# - Unix only
if(UNIX)
  # - Motif UI/Vis
  # May be required by Coin Vis driver
  # Selection also enables ToolsSG driver Xt backend
  option(GEANT4_USE_XM "Build Geant4 with Motif (X11) support" OFF)
  if((GEANT4_USE_INVENTOR AND NOT GEANT4_USE_INVENTOR_QT)
      AND NOT GEANT4_USE_XM)
    set(GEANT4_USE_XM ON CACHE BOOL "Build Geant4 with Motif (X11) support" FORCE)
    message(STATUS "Forcing GEANT4_USE_XM to ON, required by Inventor driver")
  endif()
  set(GEANT4_USE_TOOLSSG_XT_GLES ${GEANT4_USE_XM})
  geant4_add_feature(GEANT4_USE_XM "Build Geant4 with Xm Support")

  # - OpenGL/X11 Vis Driver
  # Selection also enables ToolsSG driver X11 backend
  option(GEANT4_USE_OPENGL_X11 "Build Geant4 OpenGL driver with X11 support" OFF)
  set(GEANT4_USE_TOOLSSG_X11_GLES ${GEANT4_USE_OPENGL_X11})
  geant4_add_feature(GEANT4_USE_OPENGL_X11 "Build Geant4 OpenGL driver with X11 support")

  # RayTracer driver with X11 support
  option(GEANT4_USE_RAYTRACER_X11 "Build RayTracer driver with X11 support" OFF)
  geant4_add_feature(GEANT4_USE_RAYTRACER_X11 "Build RayTracer driver with X11 support")
endif()

# Windows only
if(WIN32)
  option(GEANT4_USE_OPENGL_WIN32 "Build OpenGL driver with Win32 support" OFF)
  set(GEANT4_USE_TOOLSSG_WINDOWS_GLES ${GEANT4_USE_OPENGL_WIN32})
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

# - Qt5/6
if(GEANT4_USE_QT)
  # Use versioned targets to support Qt5 < 5.15
  # Once 5.15 is the minimum version, the "Qt${QT_VERSION_MAJOR}_..." variables can be dropped
  # - https://doc.qt.io/qt-6/cmake-manual.html
  find_package(QT NAMES Qt5 COMPONENTS Core REQUIRED)
  find_package(Qt${QT_VERSION_MAJOR} COMPONENTS Core Gui Widgets OpenGL PrintSupport REQUIRED)

  geant4_save_package_variables(Qt${QT_VERSION_MAJOR}
    Qt${QT_VERSION_MAJOR}Core_DIR
    Qt${QT_VERSION_MAJOR}Gui_DIR
    Qt${QT_VERSION_MAJOR}Widgets_DIR
    Qt${QT_VERSION_MAJOR}OpenGL_DIR
    Qt${QT_VERSION_MAJOR}PrintSupport_DIR)

  get_target_property(QT_QMAKE_EXECUTABLE ${Qt${QT_VERSION_MAJOR}Core_QMAKE_EXECUTABLE} IMPORTED_LOCATION)
  geant4_add_feature(GEANT4_USE_QT "Build Geant4 with Qt${QT_VERSION_MAJOR} support")

  # Qt3D is only supported on 5.15 and above, but always on if available
  set(QT3D_MINIMUM_VERSION 5.15.0)
  set(GEANT4_USE_QT3D OFF)

  if(QT_VERSION VERSION_GREATER_EQUAL QT3D_MINIMUM_VERSION)
    find_package(Qt${QT_VERSION_MAJOR}3DCore ${QT_VERSION} EXACT QUIET)
    find_package(Qt${QT_VERSION_MAJOR}3DExtras ${QT_VERSION} EXACT QUIET)
    find_package(Qt${QT_VERSION_MAJOR}3DRender ${QT_VERSION} EXACT QUIET)

    # Forward correct minimum version to CMake/etc files
    set(QT3D_MINIMUM_VERSION ${QT_VERSION})

    if(Qt${QT_VERSION_MAJOR}3DCore_FOUND AND Qt${QT_VERSION_MAJOR}3DExtras_FOUND AND Qt${QT_VERSION_MAJOR}3DRender_FOUND)
      set(GEANT4_USE_QT3D ON)
      geant4_save_package_variables(Qt${QT_VERSION_MAJOR}
        Qt${QT_VERSION_MAJOR}3DCore_DIR
        Qt${QT_VERSION_MAJOR}3DExtras_DIR
        Qt${QT_VERSION_MAJOR}3DRender_DIR)
      geant4_add_feature(GEANT4_USE_QT3D "Build Geant4 Qt3D driver")
    else()
      set(_g4_qt3d_missing)
      if(NOT Qt${QT_VERSION_MAJOR}3DCore_FOUND)
        list(APPEND _g4_qt3d_missing "Qt${QT_VERSION_MAJOR}3DCore")
      endif()
      if(NOT Qt${QT_VERSION_MAJOR}3DExtras_FOUND)
        list(APPEND _g4_qt3d_missing "Qt${QT_VERSION_MAJOR}3DExtras")
      endif()
      if(NOT Qt${QT_VERSION_MAJOR}3DRender_FOUND)
        list(APPEND _g4_qt3d_missing "Qt${QT_VERSION_MAJOR}3DRender")
      endif()

      message(STATUS "Disabling Geant4 Qt3D driver, missing Qt packages: ${_g4_qt3d_missing}")
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
   OR GEANT4_USE_OPENGL_WIN32)
  find_package(OpenGL REQUIRED)
  geant4_save_package_variables(OpenGL OPENGL_INCLUDE_DIR OPENGL_gl_LIBRARY)

  # X11 drivers on macOS need XQuartzGL
  if(APPLE AND
     (GEANT4_USE_INVENTOR_XT
      OR GEANT4_USE_OPENGL_X11
      OR GEANT4_USE_XM))
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
     OR GEANT4_USE_INVENTOR_XT)
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
