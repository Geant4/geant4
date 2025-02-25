#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(cellularPhantom)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Dowload geometry data file
set(GEOMETRY_NEEDS_DOWNLOAD TRUE)
set(GEOMETRY_NEEDS_UNPACK_DELETE TRUE)
set(GEOMETRY_FILE_NAME "phantoms.tar.gz")
set(GEOMETRY_FOlDER_NAME "phantoms")
set(GEOMETRY_LOCAL_FILENAME "${PROJECT_BINARY_DIR}/${GEOMETRY_FILE_NAME}")
set(GEOMETRY_DATASETS_URL
        "https://cern.ch/geant4-data/datasets/examples/advanced/dna/cellularPhantom/0/${GEOMETRY_FILE_NAME}")
set(HASH_MD5 "b663329eaa7d93396689506a798a4577")

if (EXISTS "${GEOMETRY_FOlDER_NAME}")
    set(GEOMETRY_NEEDS_DOWNLOAD FALSE)
endif ()


if (GEOMETRY_NEEDS_DOWNLOAD)
    message(STATUS "phantoms-data: attempting download: ${GEOMETRY_DATASETS_URL} ...")
    file(DOWNLOAD "${GEOMETRY_DATASETS_URL}" "${GEOMETRY_LOCAL_FILENAME}"
            INACTIVITY_TIMEOUT 500
            TIMEOUT 500
            STATUS DownloadStatus
            )

    list(GET DownloadStatus 0 DownloadReturnStatus)
    if (DownloadReturnStatus)
        message(FATAL_ERROR "phantoms-data: download FAILED: ${DownloadReturnStatus},
    This example needs internet for the phantoms data file,
    even configuring done and complied.
    Please, check your connection.
    ")
    else ()
        message(STATUS "phantoms-data: download OK")
    endif ()
endif ()


if (EXISTS "${GEOMETRY_FOlDER_NAME}")
    set(GEOMETRY_NEEDS_UNPACK_DELETE FALSE)
endif ()

if (GEOMETRY_NEEDS_UNPACK_DELETE)
    message(STATUS "Going to unpack: phantoms.tar.gz")
    execute_process(
            COMMAND ${CMAKE_COMMAND} -E tar xfz "${GEOMETRY_LOCAL_FILENAME}"
            OUTPUT_QUIET
            RESULT_VARIABLE __phantoms_untar_result
    )
    if (__phantoms_untar_result)
        message(FATAL_ERROR "phantoms-data: failed to untar file : ${GEOMETRY_LOCAL_FILENAME}")
    else ()
        message(STATUS "phantoms-data: untarred in '${PROJECT_BINARY_DIR}/phantoms' OK")
    endif ()
    message(STATUS "Going to delete: ${GEOMETRY_LOCAL_FILENAME}")
    execute_process(
            COMMAND rm "${GEOMETRY_LOCAL_FILENAME}"
    )
endif ()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(cellularPhantom cellularPhantom.cc ${sources} ${headers})
target_link_libraries(cellularPhantom ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build cellule. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(cellule_SCRIPTS
    vis.mac run.mac plot.C
  )

foreach(_script ${cellule_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS cellularPhantom DESTINATION bin)
