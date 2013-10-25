#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test47)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
find_package(Geant4 REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

find_package(ROOT QUIET)
if(NOT ROOT_FOUND)
  message(STATUS "G4 TESTS: ROOT package not found. ")  
  return()
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
		    ${ROOTSYS}/include)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(test47 test47.cc ${sources} ${headers})
target_link_libraries(test47 ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} )

# copy over the ctest config
#
set(test47_SCRIPTS ctest47-itep-proton.config)
foreach(_script ${test47_SCRIPTS})
   configure_file(${PROJECT_SOURCE_DIR}/${_script} ${PROJECT_BINARY_DIR}/${_script} COPYONLY)
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS test47 DESTINATION bin)
