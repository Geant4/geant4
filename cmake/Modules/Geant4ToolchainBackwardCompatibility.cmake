# - Script for setting up backward compatibility with GNU make user toolchain
#
# The GNU make based buildsystem for Geant4 provides a toolchain for users
# building simple Geant4 applications. The old style build and install of 
# Geant4 provides a customized set of non-standard install paths with use of 
# the toolchain dependent on environment variables pointing to the install
# paths.
#
# This script processes information on the CMake install paths, system and
# compiler to determine the following variables for backward compatibility
#
#     GEANT4_SYSTEM  : Old style system name, e.g. 'Linux', 'Darwin' or 'WIN32'
#
#     GEANT4_COMPILER: Old system compiler id, e.g. 'g++', 'VC'. 
#
#     G4INSTALL      : Location of config subdirectory containing GNU make
#                      fragment files for toolchain.
#
#     G4INCLUDE      : Old style path to location of Geant4 headers
#
#     G4LIB          : Old style library directory path. Rather than
#                      containing the actual libraries, it is expected to
#                      contain subdirectories named
#                      GEANT4_SYSTEM-GEANT4_COMPILER
#
#
# These variables are used in a CMake configuration file which is used to 
# generate shell scripts (C and Bourne flavour) the user can source to set 
# up their environment for use of the old toolchain. 
# These replace the old style 'env.(c)sh' scripts to allow users to work with 
# the new CMake built libraries transparently if their application relies on 
# the old style toolchain.
#
# The scripts are generated for both the build and install trees so that
# developers wishing to write test applications do not have to install their
# fresh build of Geant4.
#
# Compatibility with the library path style:
#
#     <prefix>/lib/G4SYSTEM-G4COMPILER
#
# is provided by installing a directory 'geant4-<version>' in the actual
# library directory and creating a symbolic link inside this directory
# pointing up one directory level.
# This will not work on Windows however, and here users are recommended to
# use Visual Studio directly, or to use the under development CMake 
# Geant4Config system.
#
# $Id: Geant4ToolchainBackwardCompatibility.cmake,v 1.10 2010-12-13 19:03:34 bmorgan Exp $
# GEANT4 Tag $Name: not supported by cvs2svn $
#

#----------------------------------------------------------------------------
# - Functions and Macros to help configuration of shell scripts.
#
#----------------------------------------------------------------------------
# MACRO(_g4tc_shell_setup)
# Set shell parameters such as program, family and common builtins for either
# bourne or cshel families
#
macro(_g4tc_shell_setup SHELL_FAMILY)
    if(${SHELL_FAMILY} STREQUAL "bourne")
        set(GEANT4_TC_SHELL_PROGRAM "/bin/sh")
        set(GEANT4_TC_SHELL_FAMILY "Bourne shell")
        set(GEANT4_TC_UNSET_COMMAND "unset")
        set(GEANT4_TC_SHELL_EXTENSION ".sh")
    elseif(${SHELL_FAMILY} STREQUAL "cshell")
        set(GEANT4_TC_SHELL_PROGRAM "/bin/csh")
        set(GEANT4_TC_SHELL_FAMILY "C shell")
        set(GEANT4_TC_UNSET_COMMAND "unsetenv")
        set(GEANT4_TC_SHELL_EXTENSION ".csh")
    else()
        # Nothing about other shells...
    endif()
endmacro()


#----------------------------------------------------------------------------
# FUNCTION(_g4tc_setenv_command)
# Return the string giving the shell command to set an environment variable
#
function(_g4tc_setenv_command TEMPLATE_NAME SHELL_FAMILY VARIABLE_NAME VARIABLE_VALUE)
    if(${SHELL_FAMILY} STREQUAL "bourne")
        set(${TEMPLATE_NAME} 
            "export ${VARIABLE_NAME}=${VARIABLE_VALUE}"
            PARENT_SCOPE
        )
    elseif(${SHELL_FAMILY} STREQUAL "cshell")
        set(${TEMPLATE_NAME} 
            "setenv ${VARIABLE_NAME} ${VARIABLE_VALUE}"
            PARENT_SCOPE
        )
    endif()
endfunction()


#----------------------------------------------------------------------------
# FUNCTION(_g4tc_prepend_path)
# return the string of commands needed to prepend a value to a path style
# environment variable
#
function(_g4tc_prepend_path TEMPLATE_NAME SHELL_FAMILY PATH_VARIABLE
    APPEND_VARIABLE)
    if(${SHELL_FAMILY} STREQUAL "bourne")
        # We have to make this section verbatim
        set(${TEMPLATE_NAME}
            "
if test \"x\$${PATH_VARIABLE}\" = \"x\" ; then
    export ${PATH_VARIABLE}=${APPEND_VARIABLE}
else
    export ${PATH_VARIABLE}=${APPEND_VARIABLE}:\${${PATH_VARIABLE}}
fi
"
            PARENT_SCOPE
        )
    elseif(${SHELL_FAMILY} STREQUAL "cshell")
        # Again, this is verbatim so final output is formatted correctly
        set(${TEMPLATE_NAME}
            "
if ( ! \${?${PATH_VARIABLE}} ) then
    setenv ${PATH_VARIABLE} ${APPEND_VARIABLE}
else
    setenv ${PATH_VARIABLE} ${APPEND_VARIABLE}:\${${PATH_VARIABLE}}
endif
"
            PARENT_SCOPE
        )
    endif()
endfunction()




#----------------------------------------------------------------------------
# MACRO(_g4tc_configure_tc_variables)
# Macro to perform the actual setting of the low level toolchain variables
# which need to be set in the final shell files.
# We do this in a separate macro so that we can wrap it in different ways for
# the install and build trees.
#
macro(_g4tc_configure_tc_variables SHELL_FAMILY)
    # - Setup the requested shell
    _g4tc_shell_setup(${SHELL_FAMILY})

    # - Standard Setup and Paths
    _g4tc_setenv_command(GEANT4_TC_G4SYSTEM ${SHELL_FAMILY} G4SYSTEM ${G4SYSTEM})
    _g4tc_setenv_command(GEANT4_TC_G4INSTALL ${SHELL_FAMILY} G4INSTALL ${G4INSTALL})
    _g4tc_setenv_command(GEANT4_TC_G4INCLUDE ${SHELL_FAMILY} G4INCLUDE ${G4INCLUDE})
    _g4tc_setenv_command(GEANT4_TC_G4LIB ${SHELL_FAMILY} G4LIB ${G4LIB})
 
    if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
        _g4tc_prepend_path(GEANT4_TC_G4LIB_PATH_SETUP ${SHELL_FAMILY} DYLD_LIBRARY_PATH ${G4LIB_DIR})
    else()
        _g4tc_prepend_path(GEANT4_TC_G4LIB_PATH_SETUP ${SHELL_FAMILY} LD_LIBRARY_PATH ${G4LIB_DIR})
    endif()

    # - Geant4 Library build setup
    # We prefer shared libs if these are built, otherwise fall back to static
    # On Win32, we also want DLLs?
    if(BUILD_SHARED_LIBS)
        _g4tc_setenv_command(GEANT4_TC_G4LIB_BUILD_SHARED ${SHELL_FAMILY} G4LIB_BUILD_SHARED 1)
        if(WIN32)
            _g4tc_setenv_command(GEANT4_TC_G4LIB_USE_DLL ${SHELL_FAMILY} G4LIB_USE_DLL 1)
        endif()
    else()
        _g4tc_setenv_command(GEANT4_TC_G4LIB_BUILD_STATIC ${SHELL_FAMILY} G4LIB_BUILD_STATIC 1)
    endif()


    # - Resource file paths
    _g4tc_setenv_command(GEANT4_TC_G4ORDPARAMTABLE_PATH ${SHELL_FAMILY} G4ORDPARAMTABLE ${G4ORDPARAMTABLE})


    # - CLHEP...
    if(GEANT4_USE_SYSTEM_CLHEP)
        # Have to use detected CLHEP paths to set base dir and others
        get_filename_component(_CLHEP_BASE_DIR ${CLHEP_INCLUDE_DIR} PATH)
        get_filename_component(_CLHEP_LIB_DIR ${CLHEP_LIBRARY} PATH)

        set(GEANT4_TC_G4LIB_USE_CLHEP "# USING SYSTEM CLHEP")
        _g4tc_setenv_command(GEANT4_TC_CLHEP_BASE_DIR ${SHELL_FAMILY} CLHEP_BASE_DIR ${_CLHEP_BASE_DIR})

        _g4tc_setenv_command(GEANT4_TC_CLHEP_INCLUDE_DIR ${SHELL_FAMILY} CLHEP_INCLUDE_DIR ${CLHEP_INCLUDE_DIR})

        _g4tc_setenv_command(GEANT4_TC_CLHEP_LIB_DIR ${SHELL_FAMILY} CLHEP_LIB_DIR ${_CLHEP_LIB_DIR})

        if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
            _g4tc_prepend_path(GEANT4_TC_CLHEP_LIB_PATH_SETUP ${SHELL_FAMILY} DYLD_LIBRARY_PATH \${CLHEP_LIB_DIR})
        else()
            _g4tc_prepend_path(GEANT4_TC_CLHEP_LIB_PATH_SETUP ${SHELL_FAMILY} LD_LIBRARY_PATH \${CLHEP_LIB_DIR})
        endif()

    else()
        # We have to configure things to point to the internal CLHEP...
        # Probably sufficient to do nothing...
        set(GEANT4_TC_G4LIB_USE_CLHEP "# USING INTERNAL CLHEP")
    endif()

    # - ZLIB...
    if(GEANT4_USE_SYSTEM_ZLIB)
        set(GEANT4_TC_G4LIB_USE_ZLIB "# USING SYSTEM ZLIB")
    else()
        _g4tc_setenv_command(GEANT4_TC_G4LIB_USE_ZLIB ${SHELL_FAMILY} G4LIB_USE_ZLIB 1)
    endif()

    # - GDML...
    if(GEANT4_USE_GDML)
        set(GEANT4_TC_G4LIB_USE_GDML "# BUILT WITH GDML SUPPORT")
    else()
        set(GEANT4_TC_G4LIB_USE_GDML "# NOT BUILT WITH GDML SUPPORT")
    endif()

    # - G3TOG4...
    if(GEANT4_USE_G3TOG4)
        _g4tc_setenv_command(GEANT4_TC_G4LIB_USE_G3TOG4 ${SHELL_FAMILY} G4LIB_USE_G3TOG4 1)
    else()
        set(GEANT4_TC_G4LIB_USE_G3TOG4 "# NOT BUILT WITH G3TOG4 SUPPORT")
    endif()

    # - USER INTERFACE AND VISUALIZATION MODULES...
    # - Terminals
    if(NOT WIN32)
        _g4tc_setenv_command(GEANT4_TC_G4UI_USE_TCSH ${SHELL_FAMILY} G4UI_USE_TCSH 1)
        set(GEANT4_TC_G4UI_USE_WIN32 "# WIN32 TERMINAL UI NOT AVAILABLE ON ${CMAKE_SYSTEM_NAME}")
    else()
        set(GEANT4_TC_G4UI_USE_TCSH "# TCSH TERMINAL UI NOT AVAILABLE ON ${CMAKE_SYSTEM_NAME}")
        _g4tc_setenv_command(GEANT4_TC_G4UI_USE_WIN32 ${SHELL_FAMILY} G4UI_USE_WIN32 1)
    endif()

    # - Qt UI AND VIS
    if(GEANT4_USE_QT)
        _g4tc_setenv_command(GEANT4_TC_G4UI_USE_QT ${SHELL_FAMILY} G4UI_USE_QT 1)
        _g4tc_setenv_command(GEANT4_TC_G4VIS_USE_OPENGLQT ${SHELL_FAMILY} G4VIS_USE_OPENGLQT 1)

        # Might need library setup, but for now recommend system install....
    else()
        set(GEANT4_TC_G4UI_USE_QT "# NOT BUILT WITH QT INTERFACE")
    endif()

    # - XM UI AND VIS
    if(GEANT4_USE_XM)
        _g4tc_setenv_command(GEANT4_TC_G4UI_USE_XM ${SHELL_FAMILY} G4UI_USE_XM 1)
        _g4tc_setenv_command(GEANT4_TC_G4VIS_USE_OPENGLXM ${SHELL_FAMILY} G4VIS_USE_OPENGLXM 1)

        # Might need library setup, but for now recommend system install....
    else()
        set(GEANT4_TC_G4UI_USE_XM "# NOT BUILT WITH XM INTERFACE")
    endif()

    # - Network DAWN
    if(GEANT4_USE_NETWORKDAWN)
        _g4tc_setenv_command(GEANT4_TC_G4VIS_USE_DAWN ${SHELL_FAMILY} G4VIS_USE_DAWN 1)
    else()
        set(GEANT4_TC_G4VIS_USE_DAWN "# NOT BUILT WITH NETWORK DAWN SUPPORT")
    endif()

    # - Network VRML
    if(GEANT4_USE_NETWORKVRML)
        _g4tc_setenv_command(GEANT4_TC_G4VIS_USE_VRML ${SHELL_FAMILY} G4VIS_USE_VRML 1)
    else()
        set(GEANT4_TC_G4VIS_USE_VRML "# NOT BUILT WITH NETWORK VRML SUPPORT")
    endif()

    # - X11 OpenGL
    if(GEANT4_USE_OPENGL_X11)
        _g4tc_setenv_command(GEANT4_TC_G4VIS_USE_OPENGLX ${SHELL_FAMILY} G4VIS_USE_OPENGLX 1)
    else()
        set(GEANT4_TC_G4VIS_USE_OPENGLX "# NOT BUILT WITH OPENGL(X11) SUPPORT")
    endif()

    # - WIN32 OpenGL
    if(GEANT4_USE_OPENGL_WIN32)
        _g4tc_setenv_command(GEANT4_TC_G4VIS_USE_OPENWIN32 ${SHELL_FAMILY} G4VIS_USE_OPENWIN32 1)
    else()
        set(GEANT4_TC_G4VIS_USE_OPENGLWIN32 "# NOT BUILT WITH OPENGL(WIN32) SUPPORT")
    endif()

    # - X11 RayTracer
    if(GEANT4_USE_RAYTRACER_X11)
        _g4tc_setenv_command(GEANT4_TC_G4VIS_USE_RAYTRACERX ${SHELL_FAMILY} G4VIS_USE_RAYTRACERX 1)
    else()
        set(GEANT4_TC_G4VIS_USE_RAYTRACERX "# NOT BUILT WITH RAYTRACER(X11) SUPPORT")
    endif()
endmacro()


#----------------------------------------------------------------------------
# MACRO(_g4tc_configure_build_tree_scripts)
# Macro to configure toolchain compatibility scripts for the build tree
#
macro(_g4tc_configure_build_tree_scripts)
    # Need to process for bourne and cshell families
    foreach(_shell bourne;cshell)
        # Generate the variables
        _g4tc_configure_tc_variables(${_shell})

        # Configure the file - goes straight into the binary dir
        configure_file(
            ${CMAKE_SOURCE_DIR}/cmake/Templates/geant4-environment-skeleton.in
            ${PROJECT_BINARY_DIR}/geant4-environment-setup${GEANT4_TC_SHELL_EXTENSION}
            @ONLY
        )
    endforeach()
endmacro()


#----------------------------------------------------------------------------
# MACRO(_g4tc_configure_install_tree_script)
# Macro to configure toolchain compatibility scripts for the install tree
#
macro(_g4tc_configure_install_tree_scripts CONFIGURE_DESTINATION INSTALL_DESTINATION)
    # Need to process for bourne and cshell families
    foreach(_shell bourne;cshell)
        # Generate the variables
        _g4tc_configure_tc_variables(${_shell})

        # Configure the file
        configure_file(
            ${CMAKE_SOURCE_DIR}/cmake/Templates/geant4-environment-skeleton.in
            ${CONFIGURE_DESTINATION}/geant4-environment-setup${GEANT4_TC_SHELL_EXTENSION}
            @ONLY
        )

        # Install it to the required location
        install(FILES
            ${CONFIGURE_DESTINATION}/geant4-environment-setup${GEANT4_TC_SHELL_EXTENSION}
            DESTINATION ${INSTALL_DESTINATION}
            PERMISSIONS 
                OWNER_READ OWNER_WRITE OWNER_EXECUTE
                GROUP_READ GROUP_EXECUTE
                WORLD_READ WORLD_EXECUTE
            COMPONENT Development
        )
    endforeach()
endmacro()


#----------------------------------------------------------------------------
# Implementation section
# 
#----------------------------------------------------------------------------
# Configure shell scripts for BUILD TREE
# This means we have to point to libraries in the build tree, but includes and
# resource files in the source tree
# This scripts never need to be relocatable
# N.B. IT WILL NOT WORK when building with VS/Xcode or any multiconfig
# buildtool because we cannot reconcile the output paths these use with those
# expected by the old toolchain...
#
set(G4SYSTEM  "${GEANT4_SYSTEM}-${GEANT4_COMPILER}")
set(G4INSTALL ${PROJECT_SOURCE_DIR})
set(G4INCLUDE ${PROJECT_SOURCE_DIR}/no_include)
set(G4LIB ${PROJECT_BINARY_DIR}/outputs/library)
set(G4LIB_DIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set(G4ORDPARAMTABLE ${G4INSTALL}/source/physics_lists/builders/OrderingParameterTable)

#----------------------------------------------------------------------------
# Configure the shell scripts for the BUILD TREE
_g4tc_configure_build_tree_scripts()


#----------------------------------------------------------------------------
# Configure shell scripts for INSTALL TREE
# This means we have to point things to their final location when installed
# These paths are all determined by the CMAKE_INSTALL_FULL directories and
# others.
# We make no attempt yet to try and make the shell scripts relocatable.
# This is likely possible, but appears to need some care.
#
#------------------------------------------------------------------------------
# Construct universal backward compatible INSTALL TREE PATHS.
#
set(G4SYSTEM  "${GEANT4_SYSTEM}-${GEANT4_COMPILER}")
set(G4INSTALL ${CMAKE_INSTALL_FULL_DATAROOTDIR}/Geant4-${Geant4_VERSION})
set(G4INCLUDE ${CMAKE_INSTALL_FULL_INCLUDEDIR}/${PROJECT_NAME})
set(G4LIB     ${CMAKE_INSTALL_FULL_LIBDIR}/Geant4-${Geant4_VERSION})
set(G4LIB_DIR ${CMAKE_INSTALL_FULL_LIBDIR})
set(G4ORDPARAMTABLE ${G4INSTALL}/physics_lists/OrderingParameterTable)


#----------------------------------------------------------------------------
# Configure the shell scripts for the INSTALL TREE
#
_g4tc_configure_install_tree_scripts(
    ${CMAKE_BINARY_DIR}/InstallTreeFiles 
    ${CMAKE_INSTALL_DATAROOTDIR}/Geant4-${Geant4_VERSION}
)


#----------------------------------------------------------------------------
# For install tree, we also need to install the config directory which contains
# all the old toolchain scripts, and to create a softlink to the G4SYSTEM
# directory.
#
install(DIRECTORY config
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}-${${PROJECT_NAME}_VERSION}
    COMPONENT Development
    FILES_MATCHING PATTERN "*.gmk"
    PATTERN "CVS" EXCLUDE
    PATTERN ".svn" EXCLUDE
    PATTERN "scripts/" EXCLUDE
)

# compatibility softlink to library directory
# NB This won't work on Windows, but should fail either.
install(CODE "execute_process(COMMAND \${CMAKE_COMMAND} -E make_directory \$ENV{DESTDIR}${CMAKE_INSTALL_FULL_LIBDIR}/Geant4-${Geant4_VERSION})")

install(CODE "execute_process(COMMAND \${CMAKE_COMMAND} -E create_symlink .. ${GEANT4_SYSTEM}-${GEANT4_COMPILER} WORKING_DIRECTORY \$ENV{DESTDIR}${CMAKE_INSTALL_FULL_LIBDIR}/Geant4-${Geant4_VERSION})")


