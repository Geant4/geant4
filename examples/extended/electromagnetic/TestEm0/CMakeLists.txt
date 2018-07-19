#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(TestEm0)

#----------------------------------------------------------------------------
# Find Geant4 package, no UI and Vis drivers activated
#
find_package(Geant4 REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})


#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executables, and link it to the Geant4 libraries
#
add_executable(TestEm0 TestEm0.cc ${sources} ${headers})
target_link_libraries(TestEm0 ${Geant4_LIBRARIES} )
add_executable(DirectAccess DirectAccess.cc ${sources} ${headers})
target_link_libraries(DirectAccess ${Geant4_LIBRARIES} )

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build TestEm0. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(TestEm0_SCRIPTS
    TestEm0.in TestEm0.out muon.mac
  )

foreach(_script ${TestEm0_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS TestEm0 DESTINATION bin)
install(TARGETS DirectAccess DESTINATION bin)

