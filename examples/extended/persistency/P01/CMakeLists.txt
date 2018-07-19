#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(P01)

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
# Find ROOT (required package)
#
find_package(ROOT REQUIRED)

# ROOT version 6 required
if(ROOT_FOUND)
  STRING(REGEX MATCH "6.*" VERSION6MATCH ${ROOT_VERSION})
  if(NOT VERSION6MATCH)
    message(FATAL_ERROR "P01 requires ROOT 6")
  endif()
endif()

# Include ROOT's CMake functions for dictionary generation
include("${ROOT_DIR}/modules/RootNewMacros.cmake")

#----------------------------------------------------------------------------
# P01 requires shared libraries
#
if(NOT Geant4_shared_FOUND)
  message(FATAL_ERROR "P01 must use shared libraries")
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIRS})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Generate dictionaries, add ROOT libraries properties
#
REFLEX_GENERATE_DICTIONARY(ExP01Classes include/ExP01Classes.hh SELECTION xml/selection.xml)
add_library(ExP01ClassesDict SHARED ${sources} ExP01Classes.cxx)
set(libsuffix .so)
set(ROOT_LIBRARY_PROPERTIES ${ROOT_LIBRARY_PROPERTIES} SUFFIX ${libsuffix})
set_target_properties(ExP01ClassesDict PROPERTIES ${ROOT_LIBRARY_PROPERTIES})
target_link_libraries(ExP01ClassesDict  ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(exampleP01 exampleP01.cc ${sources} ${headers})
target_link_libraries(exampleP01 ExP01ClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
add_executable(readHits readHits.cc  ${sources} ${headers})
target_link_libraries(readHits ExP01ClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build P01. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(P01_SCRIPTS
    run.mac vis.mac
  )

foreach(_script ${P01_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(P01 DEPENDS exampleP01 readHits )

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS exampleP01 DESTINATION bin)
install(TARGETS readHits DESTINATION bin)
install(TARGETS ExP01ClassesDict DESTINATION lib)
install(FILES ${PROJECT_BINARY_DIR}/ExP01Classes_rdict.pcm DESTINATION lib)
