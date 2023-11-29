# - Configure and install Geant4 Data Libraries
# The following functions and macros are defined to help define, install,
# reuse and export datasets and associated variables.
#
# FUNCTIONS AND MACROS
# ====================
# Public API
# ----------
# function geant4_get_datasetnames(<output variable>)
#          Store list of currently known dataset names in output variable
#
# function geant4_tupleize_datasets(<output variable>)
#          Set output variable to list of old-style dataset tuples.
#          A tuple has the format:
#            NAME/VERSION/FILENAME/EXTENSION/ENVVAR/MD5SUM
#          Provided for backward compatibility.
#
# function geant4_export_datasets([BUILD|INSTALL] <output variable>)
#          Set output variable to list of dataset tuples for export to
#          configuration scripts
#          A tuple has the format:
#            NAME/ENVVAR/PATH/FILENAME/MD5SUM
#          BUILD will set the PATH entry to the path to the dataset used
#          for the build of Geant4.
#
#          INSTALL will set the PATH entry to the path to the dataset used
#          by an install of Geant4.
#
# function geant4_add_dataset(NAME      <id>
#                             VERSION   <ver>
#                             FILENAME  <file>
#                             EXTENSION <ext>
#                             ENVVAR    <varname>
#                             MD5SUM    <md5>)
#          Define a new dataset with the following properties
#            NAME         common name of the dataset
#            VERSION      dataset version
#            FILENAME     name of packaged dataset file
#            EXTENSION    extension of packaged dataset file
#            ENVVAR       environment variable used by Geant4 code
#                         to locate this dataset
#            MD5SUM       MD5 hash of packaged dataset file
#
# function geant4_has_dataset(<name> <output variable>)
#          Set output variable to TRUE if the named dataset exists
#          in the list of known datasets.
#
# function geant4_get_dataset_property(<name> <property> <output variable>)
#          Set output variable to value of dataset "name"'s property.
#          If property does not exist, then value will be blank.
#
# function geant4_set_dataset_property(<name> <property> <value>)
#          Set value of dataset property to supplied value
#
# function geant4_configure_datasets(DESTINATION <dir>
#                                    DOWNLOAD    <installmissing>
#                                    TIMEOUT     <seconds>
#                                    )
#          Perform the actual heavy lifting to configure each declared
#          dataset for reuse or download as needed.
#
# function geant4_dataset_isinstalled(<name>
#                                     <root directory>
#                                     <output variable>)
#          Check if named dataset is installed under the root directory.
#          Set output variable to TRUE if it is, FALSE otherwise.
#
# function geant4_install_dataset(<name> <destination> <timeout>)
#          Download dataset with name to build directory, timing out the
#          download after timeout seconds, and install it
#          into its directory under the destination directory.
#
# function geant4_reuse_dataset(<name> <destination> <is installed>)
#          Reuse the dataset with name located at destination directory.
#          If it is not installed, warn user that it will need installing
#          manually in destination.
#

#-----------------------------------------------------------------------
# GEANT4 PHYSICS DATA - GLOBAL CMAKE VARIABLES
#-----------------------------------------------------------------------
# URLs, directories and dataset entries
# We may want these as properties so we can have a small API for
# retrieving them globally
#-----------------------------------------------------------------------
# Geant4 Data Repository
set(GEANT4_DATASETS_URL "https://cern.ch/geant4-data/datasets")

# Where to install data in the build tree
set(GEANT4_BUILD_FULL_DATADIR ${PROJECT_BINARY_DIR}/data)

# Where to install data in the install tree (a Default)
set(GEANT4_INSTALL_DATADIR_DEFAULT "${CMAKE_INSTALL_DATADIR}/data")

# File containing dataset list
set(GEANT4_DATASETS_DEFINITIONS "G4DatasetDefinitions")


#-----------------------------------------------------------------------
# GEANT4 PHYSICS DATA - PUBLIC CMAKE API FOR DATASET HANDLING
#-----------------------------------------------------------------------
# Properties? Shouldn't clash with the tuplized variable...
define_property(GLOBAL PROPERTY "GEANT4_DATASETS"
  BRIEF_DOCS "List of all defined Geant4 dataset names"
  FULL_DOCS
  "Each element of the list gives the name defined for the dataset.
   This name can be used in other Geant4 Data API functions to
   extract other properties of the dataset"
  )

#-----------------------------------------------------------------------
# function geant4_get_datasetnames(<output variable>)
#          Store list of currently known dataset names in output variable
#
function(geant4_get_datasetnames _output)
  get_property(_tmp GLOBAL PROPERTY GEANT4_DATASETS)
  set(${_output} ${_tmp} PARENT_SCOPE)
endfunction()

#-----------------------------------------------------------------------
# function geant4_tupleize_datasets(<output variable>)
#          Set output variable to list of old-style dataset tuples.
#          A tuple has the format:
#            NAME/VERSION/FILENAME/EXTENSION/ENVVAR/MD5SUM
#          Provided for backward compatibility.
#
function(geant4_tupleize_datasets _output)
  geant4_get_datasetnames(_names)
  set(_tmplist)

  foreach(_ds ${_names})
    set(_tuple ${_ds})
    foreach(_p VERSION FILENAME EXTENSION ENVVAR MD5SUM)
      get_property(_tmpprop GLOBAL PROPERTY ${_ds}_${_p})
      list(APPEND _tuple ${_tmpprop})
    endforeach()
    string(REPLACE ";" "/" _tuple "${_tuple}")
    list(APPEND _tmplist "${_tuple}")
  endforeach()

  set(${_output} ${_tmplist} PARENT_SCOPE)
endfunction()

#-----------------------------------------------------------------------
# function geant4_export_datasets([BUILD|INSTALL] <output variable>)
#          Set output variable to list of dataset tuples for export to
#          Geant4Config.cmake
#          A tuple has the format:
#            NAME/ENVVAR/PATH/FILENAME/MD5SUM
#          BUILD will set the PATH entry to the path to the dataset used
#          for the build of Geant4.
#
#          INSTALL will set the PATH entry to the path to the dataset used
#          by an install of Geant4.
#
#
function(geant4_export_datasets _type _output)
  geant4_get_datasetnames(_names)
  set(_tmplist)

  foreach(_ds ${_names})
    set(_tuple ${_ds})
    get_property(_tmpprop GLOBAL PROPERTY ${_ds}_ENVVAR)
    list(APPEND _tuple ${_tmpprop})

    if(${_type} STREQUAL "BUILD")
      get_property(_tmpprop GLOBAL PROPERTY ${_ds}_BUILD_DIR)
    elseif(${_type} STREQUAL "INSTALL")
      get_property(_tmpprop GLOBAL PROPERTY ${_ds}_INSTALL_DIR)
    else()
      message(FATAL_ERROR "incorrect argument to geant4_export_datasets")
    endif()
    # Ensure CMake paths
    file(TO_CMAKE_PATH "${_tmpprop}" _tmpprop)
    list(APPEND _tuple ${_tmpprop})

    get_property(_fname GLOBAL PROPERTY ${_ds}_FILENAME)
    get_property(_fvers GLOBAL PROPERTY ${_ds}_VERSION)
    get_property(_fextn GLOBAL PROPERTY ${_ds}_EXTENSION)
    list(APPEND _tuple "${_fname}.${_fvers}.${_fextn}")

    get_property(_tmpprop GLOBAL PROPERTY ${_ds}_MD5SUM)
    list(APPEND _tuple ${_tmpprop})

    # Because we have paths, use tuple separator that should not
    # appear in a path.
    string(REPLACE ";" "|" _tuple "${_tuple}")
    list(APPEND _tmplist "${_tuple}")
  endforeach()

  set(${_output} ${_tmplist} PARENT_SCOPE)
endfunction()

#-----------------------------------------------------------------------
# function geant4_add_dataset(NAME      <id>
#                             VERSION   <ver>
#                             FILENAME  <file>
#                             EXTENSION <ext>
#                             ENVVAR    <varname>
#                             MD5SUM    <md5>)
#          Define a new dataset with the following properties
#            NAME         common name of the dataset
#            VERSION      dataset version
#            FILENAME     name of packaged dataset file
#            EXTENSION    extension of packaged dataset file
#            ENVVAR       environment variable used by Geant4 code
#                         to locate this dataset
#            MD5SUM       MD5 hash of packaged dataset file
#
function(geant4_add_dataset)
  # - Parse arguments and create variables
  set(oneValueArgs NAME VERSION FILENAME EXTENSION ENVVAR MD5SUM)
  cmake_parse_arguments(_G4ADDDATA "" "${oneValueArgs}" "" ${ARGN})

  # - Fail if already defined
  geant4_has_dataset(${_G4ADDDATA_NAME} _dsexists)
  if(_dsexists)
    message(FATAL_ERROR "Dataset ${_G4ADDDATA_NAME} is already defined")
  endif()

  # - Set properties as global props <NAME>_<PROP>
  # Simple properties...
  foreach(_prop VERSION FILENAME EXTENSION ENVVAR MD5SUM)
    set_property(GLOBAL PROPERTY ${_G4ADDDATA_NAME}_${_prop} ${_G4ADDDATA_${_prop}})
  endforeach()

  # Derived properties...
  # FILE : the full filename of the packed dataset
  set_property(GLOBAL PROPERTY ${_G4ADDDATA_NAME}_FILE
    "${_G4ADDDATA_FILENAME}.${_G4ADDDATA_VERSION}.${_G4ADDDATA_EXTENSION}"
    )
  # DIRECTORY : the name of the directory that results from unpacking
  #             the packed dataset file.
  set_property(GLOBAL PROPERTY ${_G4ADDDATA_NAME}_DIRECTORY
    "${_G4ADDDATA_NAME}${_G4ADDDATA_VERSION}"
    )
  # URL : remote and full URL to the packed dataset file
  set_property(GLOBAL PROPERTY ${_G4ADDDATA_NAME}_URL
    "${GEANT4_DATASETS_URL}/${_G4ADDDATA_FILENAME}.${_G4ADDDATA_VERSION}.${_G4ADDDATA_EXTENSION}"
    )

  # - add it to the list of defined datasets
  set_property(GLOBAL APPEND PROPERTY GEANT4_DATASETS ${_G4ADDDATA_NAME})
endfunction()

#-----------------------------------------------------------------------
# function geant4_has_dataset(<name> <output variable>)
#          Set output variable to TRUE if the named dataset exists
#          in the list of known datasets.
#
function(geant4_has_dataset _name _output)
  get_property(_dslist GLOBAL PROPERTY GEANT4_DATASETS)
  list(FIND _dslist ${_name} _index)
  if(_index GREATER -1)
    set(${_output} TRUE PARENT_SCOPE)
  else()
    set(${_output} FALSE PARENT_SCOPE)
  endif()
endfunction()

#-----------------------------------------------------------------------
# function geant4_get_dataset_property(<name> <property> <output variable>)
#          Set output variable to value of dataset "name"'s property.
#          If property does not exist, then value will be blank.
#
function(geant4_get_dataset_property _name _prop _output)
  geant4_has_dataset(${_name} _dsexists)
  if(NOT _dsexists)
    message(FATAL_ERROR "non-existent dataset ${_name}")
  endif()

  get_property(_tmp GLOBAL PROPERTY ${_name}_${_prop})
  set(${_output} ${_tmp} PARENT_SCOPE)
endfunction()

#-----------------------------------------------------------------------
# function geant4_set_dataset_property(<name> <property> <value>)
#          Set value of dataset property to supplied value
#
function(geant4_set_dataset_property _name _prop _value)
  geant4_has_dataset(${_name} _dsexists)
  if(NOT _dsexists)
    message(FATAL_ERROR "non-existent dataset ${_name}")
  endif()
  set_property(GLOBAL PROPERTY ${_name}_${_prop} "${_value}")
endfunction()

#-----------------------------------------------------------------------
# function geant4_configure_datasets(DESTINATION <dir>
#                                    DOWNLOAD    <installmissing>
#                                    TIMEOUT     <seconds>
#                                    )
#          Perform the actual heavy lifting to configure each declared
#          dataset for reuse or download as needed.
#
function(geant4_configure_datasets)
  # - Parse arguments and create variables
  set(oneValueArgs DESTINATION DOWNLOAD TIMEOUT)
  cmake_parse_arguments(_G4CFGDSS "" "${oneValueArgs}" "" ${ARGN})

  # - Load configuration
  include(${GEANT4_DATASETS_DEFINITIONS})
  geant4_get_datasetnames(_dsnames)
  set(_notinstalled )

  foreach(_ds ${_dsnames})
    geant4_dataset_isinstalled(${_ds} "${_G4CFGDSS_DESTINATION}" _installed)
    if(NOT _installed AND _G4CFGDSS_DOWNLOAD)
      geant4_install_dataset(${_ds} "${_G4CFGDSS_DESTINATION}" ${_G4CFGDSS_TIMEOUT})
    else()
      geant4_reuse_dataset(${_ds} "${_G4CFGDSS_DESTINATION}" ${_installed})
      if(NOT _installed)
        list(APPEND _notinstalled ${_ds})
      endif()
    endif()
  endforeach()

  # - Produce report on datasets needing manual install, advising
  # user on how to handle these.
  # Yes, it's long, but at least it's clear :-)
  if(_notinstalled)
    message("  *WARNING*")
    message("    Geant4 has been pre-configured to look for datasets")
    message("    in the directory:")
    message(" ")
    message("    ${_G4CFGDSS_DESTINATION}")
    message(" ")
    message("    but the following datasets are NOT present on disk at")
    message("    that location:")
    message(" ")
    foreach(_missing ${_notinstalled})
      geant4_get_dataset_property(${_missing} VERSION _vers)
      message("    ${_missing} (${_vers})")
    endforeach()
    message(" ")
    message("    If you want to have these datasets installed automatically")
    message("    simply re-run cmake and set the GEANT4_INSTALL_DATA")
    message("    variable to ON. This will configure the build to download")
    message("    and install these datasets for you. For example, on the")
    message("    command line, do:")
    message(" ")
    message("    cmake -DGEANT4_INSTALL_DATA=ON <otherargs>")
    message(" ")
    message("    The variable can also be toggled in ccmake or cmake-gui.")
    message("    If you're running on a Windows system, this is the best")
    message("    solution as CMake will unpack the datasets for you")
    message("    without any further software being required")
    message(" ")
    message("    Alternatively, you can install these datasets manually")
    message("    now or after you have installed Geant4. To do this,")
    message("    download the following files:")
    message(" ")
    foreach(_missing ${_notinstalled})
      geant4_get_dataset_property(${_missing} URL _url)
      message("    ${_url}")
    endforeach()
    message(" ")
    message("    and unpack them under the directory:")
    message(" ")
    message("    ${_G4CFGDSS_DESTINATION}")
    message(" ")
    message("    As we supply the datasets packed in gzipped tar files,")
    message("    you will need the 'tar' utility to unpack them.")
    message(" ")
    message("    Nota bene: Missing datasets will not affect or break")
    message("               compilation and installation of the Geant4")
    message("               libraries.")
    message(" ")
  endif()
endfunction()

#-----------------------------------------------------------------------
# function geant4_dataset_isinstalled(<name>
#                                     <root directory>
#                                     <output variable>)
#          Check if named dataset is installed under the root directory.
#          Set output variable to TRUE if it is, FALSE otherwise.
#
function(geant4_dataset_isinstalled _name _rdirectory _output)
  geant4_get_dataset_property(${_name} DIRECTORY _dsdir)
  set(_expectedpath ${_rdirectory}/${_dsdir})

  if(IS_DIRECTORY ${_expectedpath})
    set(${_output} TRUE PARENT_SCOPE)
  else()
    set(${_output} FALSE PARENT_SCOPE)
  endif()
endfunction()

#-----------------------------------------------------------------------
# function geant4_install_dataset(<name> <destination> <timeout>)
#          Download dataset with name to build directory, timing out the
#          download after timeout seconds, and install it
#          into its directory under the destination directory.
#
function(geant4_install_dataset _name _destination _timeout)
  # - Extract needed dataset properties
  geant4_get_dataset_property(${_name} DIRECTORY _ds_dir)
  geant4_get_dataset_property(${_name} VERSION _ds_version)
  geant4_get_dataset_property(${_name} URL _ds_url)
  geant4_get_dataset_property(${_name} MD5SUM _ds_md5sum)

  message(STATUS "Configuring download of missing dataset ${_name} (${_ds_version})")

  # - Dispatch to ExternalProject 
  include(ExternalProject)

  # - We want to retain timestamps in extracted archives as the data files
  #   are "constants". CMake >= 3.24 changed this behaviour, so we setup
  #   ExternalProject based on running version
  set(__ep_extract_timestamp_arg )
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
    set(__ep_extract_timestamp_arg DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
  endif()

  ExternalProject_add(${_name}
    PREFIX Externals/${_name}-${_ds_version}
    SOURCE_DIR ${GEANT4_BUILD_FULL_DATADIR}/${_ds_dir}
    URL ${_ds_url}
    URL_MD5 ${_ds_md5sum}
    TIMEOUT ${_timeout}
    ${__ep_extract_timestamp_arg}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

  # - Configure the dataset's build and install locations
  geant4_set_dataset_property(${_name} BUILD_DIR "${PROJECT_BINARY_DIR}/data/${_ds_dir}")
  geant4_set_dataset_property(${_name} INSTALL_DIR "${_destination}/${_ds_dir}")

  # - Add install target, and report back paths...
  install(DIRECTORY ${PROJECT_BINARY_DIR}/data/${_ds_dir}
    DESTINATION ${_destination}
    COMPONENT Data
    )
endfunction()

#-----------------------------------------------------------------------
# function geant4_reuse_dataset(<name> <destination> <is installed>)
#          Reuse the dataset with name located at destination directory.
#          If it is not installed, warn user that it will need installing
#          manually in destination.
#
function(geant4_reuse_dataset _name _destination _ispresent)
  geant4_get_dataset_property(${_name} VERSION _ds_ver)
  geant4_get_dataset_property(${_name} DIRECTORY _ds_dir)
  if(_ispresent)
    message(STATUS "Reusing dataset ${_name} (${_ds_ver})")
  else()
    message(STATUS "Pre-configuring dataset ${_name} (${_ds_ver})")
  endif()

  # - In both cases, the build and install dirs are identical
  geant4_set_dataset_property(${_name} BUILD_DIR "${_destination}/${_ds_dir}")
  geant4_set_dataset_property(${_name} INSTALL_DIR "${_destination}/${_ds_dir}")
endfunction()

#-----------------------------------------------------------------------
# GEANT4 PHYSICS DATA - USER INTERFACE AND PROCESSING
#-----------------------------------------------------------------------
# User options for installing data
# - Choose a directory under which to install the data.
# - Choose whether to download and install missing datasets.
# - Change download timeout for problematic connections
#-----------------------------------------------------------------------
# Choose Physics Data Install Dir
# This follows the pattern for interface and setting as in GNUInstallDirs
if(NOT GEANT4_INSTALL_DATADIR)
  set(GEANT4_INSTALL_DATADIR "" CACHE PATH "read-only architecture independent Geant4 physics data (DATADIR/data")
  set(GEANT4_INSTALL_DATADIR "${GEANT4_INSTALL_DATADIR_DEFAULT}")
endif()

if(NOT IS_ABSOLUTE ${GEANT4_INSTALL_DATADIR})
  set(GEANT4_INSTALL_FULL_DATADIR "${CMAKE_INSTALL_PREFIX}/${GEANT4_INSTALL_DATADIR}")
else()
  set(GEANT4_INSTALL_FULL_DATADIR "${GEANT4_INSTALL_DATADIR}")
endif()

#-----------------------------------------------------------------------
# Select whether to download and install missing datasets
option(GEANT4_INSTALL_DATA "Download/Install datasets missing from GEANT4_INSTALL_DATADIR" OFF)

#-----------------------------------------------------------------------
# Provide an option for increasing the download timeout
# Helps with large datasets over slow connections.
set(GEANT4_INSTALL_DATA_TIMEOUT 1500 CACHE STRING "Timeout for Data Library downloads")
mark_as_advanced(GEANT4_INSTALL_DATA_TIMEOUT)

#-----------------------------------------------------------------------
# Set up check, download and install of needed data
#
geant4_configure_datasets(
  DESTINATION ${GEANT4_INSTALL_FULL_DATADIR}
  DOWNLOAD    ${GEANT4_INSTALL_DATA}
  TIMEOUT     ${GEANT4_INSTALL_DATA_TIMEOUT}
  )

