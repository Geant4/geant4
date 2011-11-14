# - Configure build options for Geant4 UI/Vis requiring Third Party Libraries
# All core Interfaces not requiring third party support are built automatically
# with platform differences taken into account. These are:
#
# UI Category : Always built on supported platforms
#  G4UIterminal : UNIX, WIN32
#    + - G4UItcsh : UNIX only.
#    + - G4UIcsh  : UNIX, WIN32
#  G4UIWin32    : WIN32
#
# UI category : Built on supported platforms when Third Party library available
#  G4UIQt  : UNIX, WIN32 (Requires : Qt4)
#  G4UIXm  : UNIX (Requires : Qt4)
#  G4UIXaw : DEPRECATED
#
# CONFLICT : On Mac OS X, can only have X11 support or Qt support, not both.
#            This is more for OpenGL, but for clarity, we simplify the
#            interfaces.
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
# Vis Category : Built on supported platforms when Third Party Library available
#  OpenGL : UNIX, WIN32...
#  OpenInventor :
#
# Rather than separate the UI and visualization options, we simply ask if
# the user wants to build support for a particular tool, e.g. Qt.
#
# Specific UI/Vis options are also handled here.
#
# These options are highly platform dependent, and there are even significant
# variations on each platform, so we split things up by platform.
# This gives the logic:
#
#  Linux:
#    X11, Xm, Qt, 
#
# $Id: Geant4InterfaceOptions.cmake,v 1.5 2010-11-26 18:52:16 bmorgan Exp $
# GEANT4 Tag $Name: not supported by cvs2svn $
#

#----------------------------------------------------------------------------
# MACRO _clear_opengl_variables
# - Clear all cached variables used in finding OpenGL to help in switching
# between Framework and X11 OpenGL.
macro(_clear_opengl_variables)
    set(OPENGL_INCLUDE_DIR OPENGL_INCLUDE_DIR-NOTFOUND)
    set(OPENGL_gl_LIBRARY OPENGL_gl_LIBRARY-NOTFOUND)
    set(OPENGL_glu_LIBRARY OPENGL_glu_LIBRARY-NOTFOUND)
endmacro()


#----------------------------------------------------------------------------
# Configure OpenInventor support (CROSSPLATFORM??)
#
option(GEANT4_USE_INVENTOR "Build Geant4 OpenInventor Visualization Driver" OFF)

if(GEANT4_USE_INVENTOR)
  find_package(Inventor REQUIRED)

  # On UNIX, we also require Xm and X11...
  if(UNIX)
    find_package(Motif REQUIRED)
    find_package(X11 REQUIRED)
  endif()

  GEANT4_ADD_FEATURE(
    GEANT4_USE_INVENTOR "Build OpenInventor Driver"
    )
endif()


#----------------------------------------------------------------------------
# Configure Qt Support if needed (CROSSPLATFORM).
# We do this first as the preferred interface, and also because on Mac OS X
# the user has to make a choice between Qt and X11 based interfaces.
# If the user chooses to use Qt, set a variable to prevent querying for
# X11 based OpenGL drivers.
#
if(UNIX AND APPLE)
    # Option for Qt first, but we're going to default to X11 because this
    # is guaranteed to be installed
    option(GEANT4_USE_QT "Build Geant4 with Qt support (DISABLES X11 support)" OFF)

    # Now, the other X11 based G/UIs are dependent on GEANT4_USE_QT being OFF
    CMAKE_DEPENDENT_OPTION(
        GEANT4_USE_OPENGL_X11 
        "Build Geant4 OpenGL driver with X11 support" OFF 
        "NOT GEANT4_USE_QT" OFF
    )

    CMAKE_DEPENDENT_OPTION(
        GEANT4_USE_XM
        "Build Geant4 with Xm (X11-Motif) support" OFF 
        "NOT GEANT4_USE_QT" OFF
    )

else()
    # On Windows and non-Apple Unices, we can always enable Qt support.
    option(GEANT4_USE_QT "Build Geant4 with Qt support" OFF)
    # On non-Apple Unices, we can also always choose X11/Xm
    if(UNIX)
        option(
            GEANT4_USE_OPENGL_X11 
            "Build Geant4 OpenGL driver with X11 support" 
            OFF
        )
        option(
            GEANT4_USE_XM
            "Build Geant4 with Xm (X11-Motif) support"
        )
    endif()
endif()


#----------------------------------------------------------------------------
# If Qt support requested, find Qt4
#
if(GEANT4_USE_QT)
    # Find and configure Qt and OpenGL - require 4
    # This is fine on Mac OS X because Qt will use Framework GL.
    # On WIN32 only, set QT_USE_IMPORTED_TARGETS. The Qt4 module will otherwise
    # set QT_LIBRARIES using the 'optimized A; debug Ad' technique. CMake
    # will then complain because when building DLLs we reset 
    # LINK_INTERFACE_LIBRARIES and this cannot be passed a link of this form.
    # This means we have to recreate the imported targets later though...
    if(WIN32)
        set(QT_USE_IMPORTED_TARGETS ON)
    endif()

    find_package(Qt4 REQUIRED COMPONENTS QtCore QtGui QtOpenGL)

    # Only clear on non WIN32 - doesn't seem to like this reset?
    if(NOT WIN32)
        _clear_opengl_variables()
    endif()
    find_package(OpenGL REQUIRED)

    # OpenGL part of Qt is in OpenGL component so mark the need to
    # add OpenGL.
    set(GEANT4_USE_OPENGL TRUE)
    GEANT4_ADD_FEATURE(GEANT4_USE_QT "Build Geant4 with Qt support")
endif()


#----------------------------------------------------------------------------
# If X11/Xm were requested, find X11, Xm and OpenGL as needed.
#
# - Find X11 and OpenGL
if(GEANT4_USE_OPENGL_X11 OR GEANT4_USE_XM)
    # Always need X11
    find_package(X11 REQUIRED)

    #--------------------------------------------------------------------
    # We also need Xmu... Just add it in here, copying pattern from 
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

    # Now for OpenGL - On non-Apple UNIX, just use find_package, otherwise
    # use a direct search
    if(UNIX AND NOT APPLE)
        find_package(OpenGL REQUIRED)
    else()
        # - In case we've previously run with Qt enabled, clear the cache
        # of OpenGL variables
        _clear_opengl_variables()

        # - This is for X11 drivers, so we DON't want Framework!
        set(CMAKE_FIND_FRAMEWORK_SAVE ${CMAKE_FIND_FRAMEWORK})
        set(CMAKE_FIND_FRAMEWORK NEVER)

        find_path(OPENGL_INCLUDE_DIR GL/gl.h
            PATHS /usr/X11R6/include
            NO_DEFAULT_PATH
            )

        find_library(OPENGL_gl_LIBRARY GL
            PATHS /usr/X11R6/lib
            NO_DEFAULT_PATH
            )

        find_library(OPENGL_glu_LIBRARY GLU
            PATHS /usr/X11R6/lib
            NO_DEFAULT_PATH
            )

        mark_as_advanced(OPENGL_INCLUDE_DIR OPENGL_gl_LIBRARY OPENGL_glu_LIBRARY)
        if(NOT OPENGL_INCLUDE_DIR)
            message(FATAL_ERROR "OpenGL GL/gl.h NOTFOUND")
        endif()

        if(NOT OPENGL_gl_LIBRARY)
            message(FATAL_ERROR "OpenGL GL Library NOTFOUND")
        endif()

        if(NOT OPENGL_glu_LIBRARY)
            message(FATAL_ERROR "OpenGL GLU Library NOTFOUND")
        endif()

        # Now set the standard variables...
        set(OPENGL_LIBRARIES ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})

        set(CMAKE_FIND_FRAMEWORK ${CMAKE_FIND_FRAMEWORK_SAVE})
    endif()

    # Both X11 and Xm are in the OpenGL component, so mark this are needed
    set(GEANT4_USE_OPENGL ON)

    # Add the OpenGL X11 driver as a feature
    GEANT4_ADD_FEATURE(
        GEANT4_USE_OPENGL_X11 "Build Geant4 OpenGL driver with X11 Support"
    )
endif()

# - Find Motif if needed
if(GEANT4_USE_XM)
    find_package(Motif REQUIRED)
    GEANT4_ADD_FEATURE(GEANT4_USE_XM "Build Geant4 with Xm Support")
endif()

#----------------------------------------------------------------------------
# Now configure other UNIX only drivers
#
if(UNIX)
    #------------------------------------------------------------------------
    # X11 Support for RayTracer driver
    # We can use this on all Unices because only X11 is used.
    option(GEANT4_USE_RAYTRACER_X11 "Build RayTracer driver with X11 support" OFF)

    if(GEANT4_USE_RAYTRACER_X11)
        # Find and configure X11...
        find_package(X11 REQUIRED)
        GEANT4_ADD_FEATURE(GEANT4_USE_RAYTRACER_X11 "Build RayTracer driver with X11 support")
    endif()

    #-------------------------------------------------------------------------
    # Support for Client/Server DAWN driver
    # mark as advanced because user should know what they're doing to use this
    option(GEANT4_USE_NETWORKDAWN "Build Dawn driver with Client/Server support" OFF)
    #
    # Possible headers checks for needed network parts?
    #
    GEANT4_ADD_FEATURE(GEANT4_USE_NETWORKDAWN "Build Dawn driver with Client/Server support")

    mark_as_advanced(GEANT4_USE_NETWORKDAWN)

    #--------------------------------------------------------------------------
    # Support for Client/Server VRML driver
    # mark as advanced because user should know what they're doing to use this
    option(GEANT4_USE_NETWORKVRML "Build VRML driver with Client/Server support" OFF)
    #
    # Possible header checks for needed network parts?
    #
    GEANT4_ADD_FEATURE(GEANT4_USE_NETWORKVRML "Build VRML driver with Client/Server support")
    mark_as_advanced(GEANT4_USE_NETWORKVRML)
endif()




#----------------------------------------------------------------------------
# Configure Windows OpenGL
#
if(WIN32)
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


