# - Setup for installing architecture independent read only files.
#
# There are only two main items to install here:
#
#  Geant4 Examples - basically everything under the 'examples' directory.
#
#  Geant4 Data Libraries - not supplied with code, these must be downloaded
#                          or installed from a local URL.
#
#-----------------------------------------------------------------------
# Install examples if requested
#
option(GEANT4_INSTALL_EXAMPLES "Install source code for Geant4 examples" OFF)
mark_as_advanced(GEANT4_INSTALL_EXAMPLES)

if(GEANT4_INSTALL_EXAMPLES)
  install(DIRECTORY examples
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Geant4-${Geant4_VERSION}
    COMPONENT Examples
    PATTERN "CVS" EXCLUDE
    PATTERN ".svn" EXCLUDE)
endif()

GEANT4_ADD_FEATURE(GEANT4_INSTALL_EXAMPLES "Will install source code for Geant4 examples")

#-----------------------------------------------------------------------
# Install data libraries if requested
#-----------------------------------------------------------------------
# GEANT4 PHYSICS DATA - API AND GLOBAL CMAKE VARIABLES
#
#-----------------------------------------------------------------------
# -- URLs, directories and dataset entries
# We may want these as properties so we can have a small API for
# retrieving them globally
#-----------------------------------------------------------------------
# Geant4 Data Repository
set(GEANT4_DATASETS_URL "http://geant4.cern.ch/support/source")

# Where to install data in the build tree
set(GEANT4_BUILD_FULL_DATADIR ${PROJECT_BINARY_DIR}/data)

# Where to install data in the install tree (a Default)
set(GEANT4_INSTALL_DATADIR_DEFAULT "${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}-${${PROJECT_NAME}_VERSION}/data")

# Define the known datasets as a list of tuples, tuple entries being
# forward slash separated. Messy.
# Tuple entries:
# 0 : Directory Name
# 1 : Version
# 2 : Filename
# 3 : Filename Extension
# 4 : Environment Variable
# 5 : Expected MD5 sum File
# 6 : (NOTIMPLEMENTEDYET) Marker for detecting existing install
#
set(GEANT4_DATASETS
  G4NDL/4.1/G4NDL/tar.gz/G4NEUTRONHPDATA/ff018eca2c2ca3bc32a096c2d72df64f
  G4EMLOW/6.31/G4EMLOW/tar.gz/G4LEDATA/4461d5d80a025df6657715bb33ede6bd
  PhotonEvaporation/2.2/G4PhotonEvaporation/tar.gz/G4LEVELGAMMADATA/8010e7ce8a92564e38dd3418e6040563
  RadioactiveDecay/3.5/G4RadioactiveDecay/tar.gz/G4RADIOACTIVEDATA/5940c239734db8edf8879ae79b26b404
  G4ABLA/3.0/G4ABLA/tar.gz/G4ABLADATA/d7049166ef74a592cb97df0ed4b757bd
  G4NEUTRONXS/1.2/G4NEUTRONXS/tar.gz/G4NEUTRONXSDATA/092634b9258c7bc387cb83557ff1df81
  G4PII/1.3/G4PII/tar.gz/G4PIIDATA/05f2471dbcdf1a2b17cbff84e8e83b37
  RealSurface/1.0/RealSurface/tar.gz/G4REALSURFACEDATA/0dde95e00fcd3bcd745804f870bb6884
  G4SAIDDATA/1.1/G4SAIDDATA/tar.gz/G4SAIDXSDATA/d88a31218fdf28455e5c5a3609f7216f
  )

#-----------------------------------------------------------------------
# -- API for downloading and installing data
#-----------------------------------------------------------------------
# function _geant4_data_project(<name>
#                               PREFIX installdir
#                               SOURCE_DIR wheretounpack
#                               URL whattodownload
#                               URL_MD5 expectedMD5ofdownload
#                               TIMEOUT timeoutafter(seconds)
function(_geant4_dataproject _name)
  # - Parse arguments and create any extra needed variables
  set(oneValueArgs PREFIX SOURCE_DIR URL URL_MD5 TIMEOUT)
  cmake_parse_arguments(_G4DATA "" "${oneValueArgs}" "" ${ARGN})
  get_filename_component(_G4DATA_FILE ${_G4DATA_URL} NAME)

  # - Write Download script
  file(WRITE ${PROJECT_BINARY_DIR}/${_G4DATA_PREFIX}/${_name}-Download.cmake "
  message(STATUS \"downloading ${_G4DATA_URL}...\")
  file(DOWNLOAD
    ${_G4DATA_URL}
    \"${_G4DATA_PREFIX}/${_G4DATA_FILE}\"
    TIMEOUT ${_G4DATA_TIMEOUT}
    STATUS _status)
    # LOG _log) - LOG with TIMEOUT fails due to curl issue, prefer TIMEOUT
  list(GET _status 0 _status_code)
  list(GET _status 1 _status_msg)
  if(NOT _status_code EQUAL 0)
    message(FATAL_ERROR \"error: downloading ${_G4DATA_URL} failed
      status_code : \${_status_code}
      status_msg  : \${_status_msg}
      #log : \${_log}\"
      )
  endif()
  ")

  # - Write Verify script
  file(WRITE ${PROJECT_BINARY_DIR}/${_G4DATA_PREFIX}/${_name}-Verify.cmake "
  message(STATUS \"verifying ${_G4DATA_FILE} ...\")
  execute_process(
    COMMAND \"${CMAKE_COMMAND}\" -E md5sum \"${PROJECT_BINARY_DIR}/${_G4DATA_PREFIX}/${_G4DATA_FILE}\"
    OUTPUT_VARIABLE ov
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _result_code
    )
  if(NOT _result_code EQUAL 0)
    message(FATAL_ERROR \"error: computing md5sum of ${_G4DATA_FILE} failed\")
  endif()
  string(REGEX MATCH \"^([0-9A-Fa-f]+)\" md5_actual \"\${ov}\")
  string(TOLOWER \"\${md5_actual}\" md5_actual)
  string(TOLOWER \"${_G4DATA_URL_MD5}\" md5)
  if(NOT \"\${md5}\" STREQUAL \"\${md5_actual}\")
    message(FATAL_ERROR \"error: md5sum of '${_G4DATA_FILE}' does not match expected value
  md5_expected: \${md5}
    md5_actual: \${md5_actual}
\")
  endif()
  message(STATUS \"verifying ${_G4DATA_FILE} ... done\")
  ")

  # - Write Unpack script
  if(${_G4DATA_FILE} MATCHES "(\\.|=)(tar\\.gz|zip)$")
  else()
    message(FATAL_ERROR "error: do not know how to extract '${_G4DATA_FILE}' -- file needs to be .tar.gz or zip")
  endif()
  file(WRITE ${PROJECT_BINARY_DIR}/${_G4DATA_PREFIX}/${_name}-Unpack.cmake "
  message(STATUS \"unpacking ${_G4DATA_FILE} ...\")
  file(MAKE_DIRECTORY ${_G4DATA_SOURCE_DIR})
  execute_process(COMMAND
    \${CMAKE_COMMAND} -E tar xfz \"${PROJECT_BINARY_DIR}/${_G4DATA_PREFIX}/${_G4DATA_FILE}\"
    WORKING_DIRECTORY ${_G4DATA_SOURCE_DIR}
    RESULT_VARIABLE rv)
  if(NOT rv EQUAL 0)
    message(STATUS \"extracting... [error clean up]\")
    file(REMOVE_RECURSE \"${_G4DATA_SOURCE_DIR}\")
    message(FATAL_ERROR \"error: extract of '${_G4DATA_FILE}' failed\")
  endif()
  message(STATUS \"unpacking ${_G4DATA_FILE} ... done\")
  ")

  # - Add custom commands for each step, each depending on the last
  foreach(_step Download Verify Unpack)
    if(_last_step)
      set(_deps DEPENDS ${PROJECT_BINARY_DIR}/${_G4DATA_PREFIX}/${_name}-${_last_step}.stamp)
    endif()

    add_custom_command(
      OUTPUT ${PROJECT_BINARY_DIR}/${_G4DATA_PREFIX}/${_name}-${_step}.stamp
      COMMENT "${_step} ${_name}"
      COMMAND "${CMAKE_COMMAND}" -P ${PROJECT_BINARY_DIR}/${_G4DATA_PREFIX}/${_name}-${_step}.cmake
      COMMAND "${CMAKE_COMMAND}" -E touch ${PROJECT_BINARY_DIR}/${_G4DATA_PREFIX}/${_name}-${_step}.stamp
      VERBATIM
      ${_deps}
      )

    set(_last_step ${_step})
  endforeach()

  # - Add the main target which will run all the above steps
  add_custom_target(${_name} ALL 
    COMMENT "Completed ${_name}"
    DEPENDS ${PROJECT_BINARY_DIR}/${_G4DATA_PREFIX}/${_name}-${_last_step}.stamp
    VERBATIM
    )
endfunction()


#-----------------------------------------------------------------------
# - Dispatch download task according to CMake version, and install data
# function geant4_do_install_data(<repourl>,
#                                 <dataset_tuple>,
#                                 <prefix>,
#                                 <download_timeout>
function(geant4_do_install_data _url _dataset _prefix _timeout)
  # Listify tuple and extract parameters
  string(REPLACE "/" ";" _tuple ${_dataset})
  list(GET _tuple 0 _name)
  list(GET _tuple 1 _version)
  list(GET _tuple 2 _filename)
  list(GET _tuple 3 _extension)
  list(GET _tuple 4 _envvar)
  list(GET _tuple 5 _md5sum)

  # - Dispatch to ExternalProject or our own implementation.
  # Use of URL_MD5 *and* TIMEOUT require CMake 2.8.2 or higher.
  if(${CMAKE_VERSION} VERSION_GREATER "2.8.1")
    include(ExternalProject)
    ExternalProject_Add(${_name}
      PREFIX Externals/${_filename}-${_version}
      SOURCE_DIR ${GEANT4_BUILD_FULL_DATADIR}/${_name}${_version}
      URL ${_url}/${_filename}.${_version}.${_extension}
      URL_MD5 ${_md5sum}
      TIMEOUT ${_timeout}
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ""
      )
  else()
    _geant4_dataproject(${_name}
      PREFIX Externals/${_filename}-${_version}
      SOURCE_DIR ${GEANT4_BUILD_FULL_DATADIR}
      URL ${_url}/${_filename}.${_version}.${_extension}
      URL_MD5 ${_md5sum}
      TIMEOUT ${_timeout}
      )
  endif()

  # - Add install target, and report back paths...
  install(DIRECTORY ${PROJECT_BINARY_DIR}/data/${_name}${_version}
    DESTINATION ${_prefix}
    COMPONENT Data
    )
endfunction()


#-----------------------------------------------------------------------
# GEANT4 PHYSICS DATA - USER INTERFACE AND PROCESSING
#
#-----------------------------------------------------------------------
# User options for installing data
# - Choose whether to install data, which both files and environment
# - Choose a directory under which to install the data.
# - Change download timeout for problematic connections
#
option(GEANT4_INSTALL_DATA "Download and install Geant4 Data Libraries" OFF)

#-----------------------------------------------------------------------
# Choose Physics Data Install Dir
# This follows the pattern for interface and setting as in GNUInstallDirs
if(NOT GEANT4_INSTALL_DATADIR)
  set(GEANT4_INSTALL_DATADIR "" CACHE PATH "read-only architecture independent Geant4 physics data (DATAROOTDIR/${GEANT4_INSTALL_DATADIR_DEFAULT}")
  set(GEANT4_INSTALL_DATADIR "${GEANT4_INSTALL_DATADIR_DEFAULT}")
endif()

if(NOT IS_ABSOLUTE ${GEANT4_INSTALL_DATADIR})
  set(GEANT4_INSTALL_FULL_DATADIR "${CMAKE_INSTALL_PREFIX}/${GEANT4_INSTALL_DATADIR}")
else()
  set(GEANT4_INSTALL_FULL_DATADIR "${GEANT4_INSTALL_DATADIR}")
endif()

mark_as_advanced(GEANT4_INSTALL_DATADIR)

#-----------------------------------------------------------------------
# Provide an option for increasing the download timeout
# Helps with large datasets over slow connections.
set(GEANT4_INSTALL_DATA_TIMEOUT 1500 CACHE STRING "Timeout for Data Library download")
mark_as_advanced(GEANT4_INSTALL_DATA_TIMEOUT)

#-----------------------------------------------------------------------
# Set up check, download and install of needed data
#
if(GEANT4_INSTALL_DATA)
  foreach(_dataset ${GEANT4_DATASETS})
    geant4_do_install_data(${GEANT4_DATASETS_URL} ${_dataset} ${GEANT4_INSTALL_FULL_DATADIR} ${GEANT4_INSTALL_DATA_TIMEOUT})
  endforeach()
  GEANT4_ADD_FEATURE(GEANT4_INSTALL_DATA "Will download and install data libraries")
endif()


