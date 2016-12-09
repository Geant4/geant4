#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(dicomReader)

#----------------------------------------------------------------------------
# Find Geant4 package, no UI and Vis drivers activated
#
find_package(Geant4 REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find DCMTK installation
#
find_package(DCMTK REQUIRED)
add_definitions(-DG4_DCMTK)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
                    ${DCMTK_INCLUDE_DIRS})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the library
#
add_library(dicomReader ${sources} ${headers})
target_link_libraries(dicomReader ${Geant4_LIBRARIES}  ${DCMTK_LIBRARIES})

#----------------------------------------------------------------------------
# Install the library to 'lib' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS dicomReader DESTINATION lib)
