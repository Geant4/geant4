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

#----------------------------------------------------------------------------
# - Functions and Macros to help configuration of shell scripts.
#
#----------------------------------------------------------------------------
# MACRO(_g4tc_shell_setup)
# Set shell parameters such as program, family and common builtins for either
# bourne or cshell families
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
# FUNCTION(_g4tc_selflocate)
# return the shell commands need to locate the directory in which the script
# is located when it is sourced.
#
function(_g4tc_selflocate TEMPLATE_NAME SHELL_FAMILY SCRIPT_NAME LOCATION_VARIABLE)
  if(${SHELL_FAMILY} STREQUAL "bourne")
    set(${TEMPLATE_NAME} 
      "# Self locate script when sourced
if [ \"x\${BASH_ARGV[0]}\" = \"x\" ]; then
  # Not bash, so rely on sourcing from correct location
  if [ ! -f ${SCRIPT_NAME}${GEANT4_TC_SHELL_EXTENSION} ]; then
    echo 'ERROR: cd to location of ${SCRIPT_NAME} script and source it there'
    return 1
  fi
  ${LOCATION_VARIABLE}=\$\(pwd\)
else
  g4sls_sourced_dir=\$\(dirname \${BASH_ARGV[0]}\)
  ${LOCATION_VARIABLE}=$\(cd \$g4sls_sourced_dir; pwd\)
fi
      "
      PARENT_SCOPE
      )

  elseif(${SHELL_FAMILY} STREQUAL "cshell")
    set(${TEMPLATE_NAME}
      "# Self locate script when sourced
set ARGS=(\$_)
set g4sls_sourced_dir=\"`dirname \${ARGS[2]}`\"
set ${LOCATION_VARIABLE}=\"`cd \${g4sls_sourced_dir}; pwd`\"
"
      PARENT_SCOPE
      )
    endif()
endfunction()




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
# FUNCTION(_g4tc_setenv_ifnotset_command)
# Return the string giving the shell comment to set an environment variable
# if the variable is not already set
#
function(_g4tc_setenv_ifnotset_command TEMPLATE_NAME SHELL_FAMILY VARIABLE_NAME VARIABLE_VALUE)
  if(${SHELL_FAMILY} STREQUAL "bourne")
    # Have to make this section verbatim to get correct formatting
    set(${TEMPLATE_NAME}
      "
if test \"x\$${VARIABLE_NAME}\" = \"x\" ; then
  export ${VARIABLE_NAME}=${VARIABLE_VALUE}
fi
"
      PARENT_SCOPE
      )
  elseif(${SHELL_FAMILY} STREQUAL "cshell")
    # Again, verbatim to get correct formatting...
    set(${TEMPLATE_NAME} 
      "
if ( ! \${?${VARIABLE_NAME}} ) then
  setenv ${VARIABLE_NAME} ${VARIABLE_VALUE}
endif
"
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
macro(_g4tc_configure_tc_variables SHELL_FAMILY SCRIPT_NAME)
  # - Setup the requested shell
  _g4tc_shell_setup(${SHELL_FAMILY})

  # - Locate self
  _g4tc_selflocate(GEANT4_TC_LOCATE_SELF_COMMAND ${SHELL_FAMILY} ${SCRIPT_NAME} geant4make_root)


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

  _g4tc_setenv_ifnotset_command(GEANT4_TC_G4WORKDIR_SETUP ${SHELL_FAMILY} G4WORKDIR ${G4WORKDIR_DEFAULT})
  _g4tc_prepend_path(GEANT4_TC_G4WORKDIR_PATH_SETUP ${SHELL_FAMILY} PATH
  \${G4WORKDIR}/bin/\${G4SYSTEM})

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
  if(GEANT4_INSTALL_DATA)
    _g4tc_setenv_command(GEANT4_TC_G4ABLADATA ${SHELL_FAMILY} G4ABLADATA ${G4ABLADATA_PATH})
    _g4tc_setenv_command(GEANT4_TC_G4LEDATA ${SHELL_FAMILY} G4LEDATA ${G4LEDATA_PATH})
    _g4tc_setenv_command(GEANT4_TC_G4LEVELGAMMADATA ${SHELL_FAMILY} G4LEVELGAMMADATA ${G4LEVELGAMMADATA_PATH})
    _g4tc_setenv_command(GEANT4_TC_G4NEUTRONHPDATA ${SHELL_FAMILY} G4NEUTRONHPDATA ${G4NEUTRONHPDATA_PATH})
    _g4tc_setenv_command(GEANT4_TC_G4NEUTRONXSDATA ${SHELL_FAMILY} G4NEUTRONXSDATA ${G4NEUTRONXSDATA_PATH})
    _g4tc_setenv_command(GEANT4_TC_G4PIIDATA ${SHELL_FAMILY} G4PIIDATA ${G4PIIDATA_PATH})
    _g4tc_setenv_command(GEANT4_TC_G4RADIOACTIVEDATA ${SHELL_FAMILY} G4RADIOACTIVEDATA ${G4RADIOACTIVEDATA_PATH})
    _g4tc_setenv_command(GEANT4_TC_G4REALSURFACEDATA ${SHELL_FAMILY} G4REALSURFACEDATA ${G4REALSURFACEDATA_PATH})
  else()
    set(GEANT4_TC_G4ABLADATA        "# ABLA Data not installed")
    set(GEANT4_TC_G4LEDATA          "# EMLOW Data not installed")
    set(GEANT4_TC_G4LEVELGAMMADATA  "# Photon Evaporation Data not installed")
    set(GEANT4_TC_G4NEUTRONHPDATA   "# NDL Data not installed")
    set(GEANT4_TC_G4NEUTRONXSDATA   "# Neutron Cross Section Data not installed")
    set(GEANT4_TC_G4PIIDATA         "# Shell Ionization Cross Section Data not installed")
    set(GEANT4_TC_G4RADIOACTIVEDATA "# Radioactive Decay Data not installed")
    set(GEANT4_TC_G4REALSURFACEDATA "# RealSurface Data not installed")
  endif()


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
    _g4tc_setenv_command(GEANT4_TC_G4LIB_USE_GDML ${SHELL_FAMILY} G4LIB_USE_GDML 1)
    # Backward compatibility requires XERCESCROOT to be set
    # As this is a 'rootdir' determine it from the XERCESC_INCLUDE_DIR
    # variable...
    get_filename_component(_xercesc_root ${XERCESC_INCLUDE_DIR} PATH)
    _g4tc_setenv_command(GEANT4_TC_GDML_PATH_SETUP ${SHELL_FAMILY} XERCESCROOT ${_xercesc_root})
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

  # - OpenInventor
  if(GEANT4_USE_INVENTOR)
    if(UNIX)
      _g4tc_setenv_command(GEANT4_TC_G4VIS_USE_OPENINVENTOR ${SHELL_FAMILY} G4VIS_USE_OIX 1)
    else()
      _g4tc_setenv_command(GEANT4_TC_G4VIS_USE_OPENINVENTOR ${SHELL_FAMILY} G4VIS_USE_OIWIN32 1)
    endif()
  else()
    set(GEANT4_TC_G4VIS_USE_OPENINVENTOR "# NOT BUILT WITH INVENTOR SUPPORT")
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
macro(_g4tc_configure_build_tree_scripts SCRIPT_NAME)
  # Need to process for bourne and cshell families
  foreach(_shell bourne;cshell)
    # Generate the variables
    _g4tc_configure_tc_variables(${_shell} ${SCRIPT_NAME})

    # Configure the file - goes straight into the binary dir
    configure_file(
      ${CMAKE_SOURCE_DIR}/cmake/Templates/geant4make-skeleton.in
      ${PROJECT_BINARY_DIR}/${SCRIPT_NAME}${GEANT4_TC_SHELL_EXTENSION}
      @ONLY
      )
  endforeach()
endmacro()


#----------------------------------------------------------------------------
# MACRO(_g4tc_configure_install_tree_script)
# Macro to configure toolchain compatibility scripts for the install tree
#
macro(_g4tc_configure_install_tree_scripts CONFIGURE_DESTINATION SCRIPT_NAME INSTALL_DESTINATION)
  # Need to process for bourne and cshell families
  foreach(_shell bourne;cshell)
    # Generate the variables
    _g4tc_configure_tc_variables(${_shell} ${SCRIPT_NAME})

    # Configure the file
    configure_file(
      ${CMAKE_SOURCE_DIR}/cmake/Templates/geant4make-skeleton.in
      ${CONFIGURE_DESTINATION}/${SCRIPT_NAME}${GEANT4_TC_SHELL_EXTENSION}
      @ONLY
      )

    # Install it to the required location
    install(FILES
      ${CONFIGURE_DESTINATION}/${SCRIPT_NAME}${GEANT4_TC_SHELL_EXTENSION}
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
# This script never needs to be relocatable, so we don't need to use the
# self location functionality.
# N.B. IT WILL NOT WORK when building with VS/Xcode or any multiconfig
# buildtool because we cannot reconcile the output paths these use with those
# expected by the old toolchain...
#
set(G4SYSTEM  "${GEANT4_SYSTEM}-${GEANT4_COMPILER}")
set(G4INSTALL ${PROJECT_SOURCE_DIR})
set(G4INCLUDE ${PROJECT_SOURCE_DIR}/no_include)
set(G4LIB ${PROJECT_BINARY_DIR}/outputs/library)
set(G4LIB_DIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set(G4WORKDIR_DEFAULT "\$HOME/geant4_workdir")

# - Data
# Messy for now...
if(GEANT4_INSTALL_DATA)
  set(_g4datarootdir ${PROJECT_BINARY_DIR}/data)

  foreach(_ds ${GEANT4_DATASETS})
    string(REPLACE "/" ";" _tuple ${_ds})
    list(GET _tuple 0 _name)
    list(GET _tuple 1 _vers)
    list(GET _tuple 4 _envvarname)

    set(${_envvarname}_PATH ${_g4datarootdir}/${_name}${_vers})
  endforeach()
endif()

#----------------------------------------------------------------------------
# Configure the shell scripts for the BUILD TREE
_g4tc_configure_build_tree_scripts(geant4make)



#----------------------------------------------------------------------------
# Configure shell scripts for INSTALL TREE
# This means we have to point things to their final location when installed
# These paths are all determined by the CMAKE_INSTALL_FULL directories and
# others.
# If we are relocatable, then the structure we will have is
# +- CMAKE_INSTALL_PREFIX
#    +- LIBDIR/Geant4-VERSION (G4LIB)
#    +- INCLUDEDIR/Geant4     (G4INCLUDE)
#    +- DATAROOTDIR/Geant4-VERSION/
#       +- geant4make              (G4INSTALL!)
#          +- geant4make.(c)sh
#          +- config/

#------------------------------------------------------------------------------
# Construct universal backward compatible INSTALL TREE PATHS.
#
set(G4SYSTEM  "${GEANT4_SYSTEM}-${GEANT4_COMPILER}")
set(G4INSTALL "\"\$geant4make_root\"")

# Now need relative paths between 'G4INSTALL' and include/lib dirs
# - Include dir
file(RELATIVE_PATH
  G4MAKE_TO_INCLUDEDIR
  ${CMAKE_INSTALL_FULL_DATAROOTDIR}/Geant4-${Geant4_VERSION}/geant4make
  ${CMAKE_INSTALL_FULL_INCLUDEDIR}/${PROJECT_NAME}
  )
set(G4INCLUDE "\"`cd \$geant4make_root/${G4MAKE_TO_INCLUDEDIR}\; pwd`\"")

# - Lib dir
file(RELATIVE_PATH
  G4MAKE_TO_LIBDIR
  ${CMAKE_INSTALL_FULL_DATAROOTDIR}/Geant4-${Geant4_VERSION}/geant4make
  ${CMAKE_INSTALL_FULL_LIBDIR}
  )
set(G4LIB "\"`cd \$geant4make_root/${G4MAKE_TO_LIBDIR}/Geant4-${Geant4_VERSION}\; pwd`\"")
set(G4LIB_DIR "\"`cd \$geant4make_root/${G4MAKE_TO_LIBDIR}\; pwd`\"")

set(G4WORKDIR_DEFAULT "\$HOME/geant4_workdir")

# - Data
# Messy for now...
if(GEANT4_INSTALL_DATA)
  file(RELATIVE_PATH
    G4MAKE_TO_DATADIR
    ${CMAKE_INSTALL_FULL_DATAROOTDIR}/Geant4-${Geant4_VERSION}/geant4make
    ${CMAKE_INSTALL_FULL_DATAROOTDIR}/Geant4-${Geant4_VERSION}/data
    )

  foreach(_ds ${GEANT4_DATASETS})
    string(REPLACE "/" ";" _tuple ${_ds})
    list(GET _tuple 0 _name)
    list(GET _tuple 1 _vers)
    list(GET _tuple 4 _envvarname)

    set(${_envvarname}_PATH "\"`cd \$geant4make_root/${G4MAKE_TO_DATADIR}/${_name}${_vers}\; pwd`\"")
  endforeach()
endif()

#----------------------------------------------------------------------------
# Configure the shell scripts for the INSTALL TREE
#
_g4tc_configure_install_tree_scripts(
    ${CMAKE_BINARY_DIR}/InstallTreeFiles
    geant4make
    ${CMAKE_INSTALL_DATAROOTDIR}/Geant4-${Geant4_VERSION}/geant4make
)


#----------------------------------------------------------------------------
# For install tree, we also need to install the config directory which contains
# all the old toolchain scripts, and to create a softlink to the G4SYSTEM
# directory.
#
install(DIRECTORY config
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Geant4-${Geant4_VERSION}/geant4make
    COMPONENT Development
    FILES_MATCHING PATTERN "*.gmk"
    PATTERN "CVS" EXCLUDE
    PATTERN ".svn" EXCLUDE
    PATTERN "scripts/" EXCLUDE
)

# compatibility softlink to library directory
# NB This won't work on Windows, but shouldn't fail either.
install(CODE "execute_process(COMMAND \${CMAKE_COMMAND} -E make_directory \$ENV{DESTDIR}${CMAKE_INSTALL_FULL_LIBDIR}/Geant4-${Geant4_VERSION})")

install(CODE "execute_process(COMMAND \${CMAKE_COMMAND} -E create_symlink .. ${GEANT4_SYSTEM}-${GEANT4_COMPILER} WORKING_DIRECTORY \$ENV{DESTDIR}${CMAKE_INSTALL_FULL_LIBDIR}/Geant4-${Geant4_VERSION})")



#----------------------------------------------------------------------------
# TEMPORARY
# Configure environment setup script for install of Geant4
# Temporarily here to keep all shell setup in one place.
# Later, should be refactored into its own module, with module containing
# all the shell tools above.
#
# - Script base name (without extension
set(_scriptbasename geant4)

# - Relative path between bindir (where script is) and library directory
file(RELATIVE_PATH
  G4ENV_BINDIR_TO_LIBDIR
  ${CMAKE_INSTALL_FULL_BINDIR}
  ${CMAKE_INSTALL_FULL_LIBDIR}
  )

# - Data
# Messy for now...
if(GEANT4_INSTALL_DATA)
  file(RELATIVE_PATH
    G4ENV_BINDIR_TO_DATADIR
    ${CMAKE_INSTALL_FULL_BINDIR}
    ${CMAKE_INSTALL_FULL_DATAROOTDIR}/Geant4-${Geant4_VERSION}/data
    )

  foreach(_ds ${GEANT4_DATASETS})
    string(REPLACE "/" ";" _tuple ${_ds})
    list(GET _tuple 0 _name)
    list(GET _tuple 1 _vers)
    list(GET _tuple 4 _envvarname)

    set(${_envvarname}_PATH "\"`cd \$geant4_envbindir/${G4ENV_BINDIR_TO_DATADIR}/${_name}${_vers}\; pwd`\"")
  endforeach()
endif()




foreach(_shell bourne;cshell)
  # Setup the shell
  _g4tc_shell_setup(${_shell})

  # Set script full name
  set(_scriptfullname ${_scriptbasename}${GEANT4_TC_SHELL_EXTENSION})

  # Set locate self command
  _g4tc_selflocate(GEANT4_ENV_SELFLOCATE_COMMAND 
    ${_shell} 
    ${_scriptbasename} 
    geant4_envbindir
    )

  # Set path, which should be where the script itself is installed
  _g4tc_prepend_path(GEANT4_ENV_BINPATH_SETUP 
    ${_shell}
    PATH
    "\"\$geant4_envbindir\""
    )

  # Set library path, based on relative paths between bindir and libdir
  if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    set(_libpathname DYLD_LIBRARY_PATH)
  else()
    set(_libpathname LD_LIBRARY_PATH)
  endif()

  _g4tc_prepend_path(GEANT4_ENV_LIBPATH_SETUP
    ${_shell}
    ${_libpathname}
    "\"`cd $geant4_envbindir/${G4ENV_BINDIR_TO_LIBDIR}; pwd`\""
    )

  # - Set data paths
  if(GEANT4_INSTALL_DATA)
    _g4tc_setenv_command(GEANT4_ENV_G4ABLADATA ${_shell} G4ABLADATA ${G4ABLADATA_PATH})
    _g4tc_setenv_command(GEANT4_ENV_G4LEDATA ${_shell} G4LEDATA ${G4LEDATA_PATH})
    _g4tc_setenv_command(GEANT4_ENV_G4LEVELGAMMADATA ${_shell} G4LEVELGAMMADATA ${G4LEVELGAMMADATA_PATH})
    _g4tc_setenv_command(GEANT4_ENV_G4NEUTRONHPDATA ${_shell} G4NEUTRONHPDATA ${G4NEUTRONHPDATA_PATH})
    _g4tc_setenv_command(GEANT4_ENV_G4NEUTRONXSDATA ${_shell} G4NEUTRONXSDATA ${G4NEUTRONXSDATA_PATH})
    _g4tc_setenv_command(GEANT4_ENV_G4PIIDATA ${_shell} G4PIIDATA ${G4PIIDATA_PATH})
    _g4tc_setenv_command(GEANT4_ENV_G4RADIOACTIVEDATA ${_shell} G4RADIOACTIVEDATA ${G4RADIOACTIVEDATA_PATH})
    _g4tc_setenv_command(GEANT4_ENV_G4REALSURFACEDATA ${_shell} G4REALSURFACEDATA ${G4REALSURFACEDATA_PATH})
  else()
    set(GEANT4_ENV_G4ABLADATA        "# ABLA Data not installed")
    set(GEANT4_ENV_G4LEDATA          "# EMLOW Data not installed")
    set(GEANT4_ENV_G4LEVELGAMMADATA  "# Photon Evaporation Data not installed")
    set(GEANT4_ENV_G4NEUTRONHPDATA   "# NDL Data not installed")
    set(GEANT4_ENV_G4NEUTRONXSDATA   "# Neutron Cross Section Data not installed")
    set(GEANT4_ENV_G4PIIDATA         "# Shell Ionization Cross Section Data not installed")
    set(GEANT4_ENV_G4RADIOACTIVEDATA "# Radioactive Decay Data not installed")
    set(GEANT4_ENV_G4REALSURFACEDATA "# RealSurface Data not installed")
  endif()


  # Configure the file
  configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/Templates/geant4-env-skeleton.in
    ${PROJECT_BINARY_DIR}/InstallTreeFiles/${_scriptfullname}
    @ONLY
    )

  # Install it to the required location
  install(FILES
    ${PROJECT_BINARY_DIR}/InstallTreeFiles/${_scriptfullname}
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    PERMISSIONS 
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE
    WORLD_READ WORLD_EXECUTE
    COMPONENT Runtime
    )
endforeach()

