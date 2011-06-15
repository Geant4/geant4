# Geant4InterfaceOptions - module to configure build options for Geant4 UI/Vis
#
# Rather than separate the UI and visualization options, we simply ask if
# the user wants to build support for a particular tool, e.g. Qt.
#
# Specific UI/Vis options are also handled here.
#
# $Id: Geant4InterfaceOptions.cmake,v 1.5 2010-11-26 18:52:16 bmorgan Exp $
# GEANT4 Tag $Name: not supported by cvs2svn $
#

if(UNIX)
    #------------------------------------------------------------------------
    # X11 (X11 OpenGL) Support (Vis Only)
    option(GEANT4_USE_OPENGL_X11 "Build Geant4 OpenGL driver with X11 support" OFF)

    if(GEANT4_USE_OPENGL_X11)
        # Find and configure X11 and OpenGL
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

        # Find OpenGL....
        find_package(OpenGL REQUIRED)

        set(GEANT4_USE_OPENGL TRUE)
        GEANT4_ADD_FEATURE(GEANT4_USE_OPENGL_X11 "Build OpenGL driver with X11 support")
    endif()

    
    #------------------------------------------------------------------------
    # XAW (X11 Athena) Support (UI only)
    #option(GEANT4_USE_XAW "Build Geant4 with XAW (X11 Athena) support" OFF)

    if(GEANT4_USE_XAW)
        # Find and configure XAW...
    endif()

    #--------------------------------------------------------------------------
    # Xm (Motif) Support (UI and Vis)
    #option(GEANT4_USE_XM "Build Geant4 with Xm (X11-Motif) support" OFF)

    if(GEANT4_USE_XM)
        # Find and configure Xm...
    endif()

    #--------------------------------------------------------------------------
    # X11 Support for RayTracer driver
    option(GEANT4_USE_RAYTRACERX "Build RayTracer driver with X11 support" OFF)

    if(GEANT4_USE_RAYTRACERX)
        # Find and configure X11...
        find_package(X11 REQUIRED)
        GEANT4_ADD_FEATURE(GEANT4_USE_RAYTRACERX "Build RayTracer driver with X11 support")
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

#------------------------------------------------------------------------------
# Qt Support (UI and Vis)
option(GEANT4_USE_QT "Build Geant4 with Qt support" OFF)

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
    find_package(OpenGL REQUIRED)
    set(GEANT4_USE_OPENGL TRUE)
    GEANT4_ADD_FEATURE(GEANT4_USE_QT "Build Geant4 with Qt support")
endif()

