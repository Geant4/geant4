#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(molecular)

#----------------------------------------------------------------------------
# Find Geant4 package
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if (WITH_GEANT4_UIVIS)
    find_package(Geant4 REQUIRED ui_all vis_all)
else ()
    find_package(Geant4 REQUIRED)
endif ()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include)
include(${Geant4_USE_FILE})


#----------------------------------------------------------------------------
# Dowload geometry data file
set(GEOMETRY_NEEDS_DOWNLOAD TRUE)
set(GEOMETRY_NEEDS_UNPACK_DELETE TRUE)
set(GEOMETRY_FILE_NAME "geometries.tar.gz")
set(GEOMETRY_FOlDER_NAME "geometries")
set(GEOMETRY_LOCAL_FILENAME "${PROJECT_BINARY_DIR}/${GEOMETRY_FILE_NAME}")
set(GEOMETRY_DATASETS_URL
        "http://geant4-dna.in2p3.fr/geometries/molecularDNA/${GEOMETRY_FILE_NAME}")
set(HASH_MD5 "3e28151dc4c4647af3ae37d0385fc443")

if (EXISTS "${GEOMETRY_FOlDER_NAME}")
    set(GEOMETRY_NEEDS_DOWNLOAD FALSE)
endif ()


if (GEOMETRY_NEEDS_DOWNLOAD)
    message(STATUS "geometries-data: attempting download: ${GEOMETRY_DATASETS_URL} ...")
    file(DOWNLOAD "${GEOMETRY_DATASETS_URL}" "${GEOMETRY_LOCAL_FILENAME}"
            INACTIVITY_TIMEOUT 500
            TIMEOUT 500
            STATUS DownloadStatus
            )

    list(GET DownloadStatus 0 DownloadReturnStatus)
    if (DownloadReturnStatus)
        message(FATAL_ERROR "geometries-data: download FAILED: ${DownloadReturnStatus},
    This example needs internet for the geometry data file,
    even configuring done and complied.
    Please, check your connection.
    ")
    else ()
        message(STATUS "geometries-data: download OK")
    endif ()
endif ()


if (EXISTS "${GEOMETRY_FOlDER_NAME}")
    set(GEOMETRY_NEEDS_UNPACK_DELETE FALSE)
endif ()

if (GEOMETRY_NEEDS_UNPACK_DELETE)
    message(STATUS "Going to unpack: geometries.tar.gz")
    execute_process(
            COMMAND ${CMAKE_COMMAND} -E tar xfz "${GEOMETRY_LOCAL_FILENAME}"
            OUTPUT_QUIET
            RESULT_VARIABLE __geometries_untar_result
    )
    if (__geometries_untar_result)
        message(FATAL_ERROR "geometries-data: failed to untar file : ${GEOMETRY_LOCAL_FILENAME}")
    else ()
        message(STATUS "geometries-data: untarred in '${PROJECT_BINARY_DIR}/geometries' OK")
    endif ()
    message(STATUS "Going to delete: ${GEOMETRY_LOCAL_FILENAME}")
    execute_process(
            COMMAND rm "${GEOMETRY_LOCAL_FILENAME}"
    )
endif ()
#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(molecular molecular.cc ${sources} ${headers})
target_link_libraries(molecular ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build molecular_proj. This is so that we can run the executable directly because 
# it relies on these scripts being in the current working directory.
#
file(GLOB molecular_SCRIPTS
        ${PROJECT_SOURCE_DIR}/*.in
        ${PROJECT_SOURCE_DIR}/*.mac
        ${PROJECT_SOURCE_DIR}/*.C)

foreach (_script ${molecular_SCRIPTS})
    configure_file(
            ${_script}
            ${PROJECT_BINARY_DIR}/.
            COPYONLY
    )
endforeach ()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS molecular DESTINATION bin)
