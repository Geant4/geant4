#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(MCTruth)

#----------------------------------------------------------------------------
# Find Geant4 package, no UI and Vis drivers activated
#
find_package(Geant4 REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})


#----------------------------------------------------------------------------
# Find HepMC (required package)
#
if (USE_HEPMC3)
 add_definitions(-DG4LIB_USE_HEPMC3)
 find_package(HepMC3 REQUIRED)
 set (EVENT_RECORD_INCLUDE_DIR ${HEPMC3_INCLUDE_DIR})
 set (EVENT_RECORD_LIBRARIES ${HEPMC3_LIB})
else()
 find_package(HepMC2 REQUIRED)
 set (EVENT_RECORD_INCLUDE_DIR ${HEPMC2_INCLUDE_DIR})
 set (EVENT_RECORD_LIBRARIES ${HEPMC2_LIBRARIES} )
endif()
#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
                    ${GENERATOR_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(mctruthex mctruthex.cc ${sources} ${headers})
target_link_libraries(mctruthex ${Geant4_LIBRARIES} ${EVENT_RECORD_LIBRARIES} )

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build MCTruth. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(MCTruth_SCRIPTS
   
  )

foreach(_script ${MCTruth_SCRIPTS})
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
add_custom_target(MCTruth DEPENDS mctruthex)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS mctruthex DESTINATION bin)

