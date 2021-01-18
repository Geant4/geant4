#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
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
# Find HepMC3 (required package)
#
add_definitions(-DG4LIB_USE_HEPMC3)
find_package(HepMC3 REQUIRED)
set (EVENT_RECORD_INCLUDE_DIR ${HEPMC3_INCLUDE_DIR})
set (EVENT_RECORD_LIBRARIES ${HEPMC3_LIB})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
                    ${GENERATOR_INCLUDE_DIR}
                    ${EVENT_RECORD_INCLUDE_DIR}
                    )
set( sources 
               ${PROJECT_SOURCE_DIR}/src/DetectorConstruction.cc         
               ${PROJECT_SOURCE_DIR}/src/MCTruthEventAction.cc           
               ${PROJECT_SOURCE_DIR}/src/MCTruthTrackInformation.cc      
               ${PROJECT_SOURCE_DIR}/src/PrimaryGeneratorAction.cc
               ${PROJECT_SOURCE_DIR}/src/MCTruthConfig.cc                
               ${PROJECT_SOURCE_DIR}/src/MCTruthManager.cc               
               ${PROJECT_SOURCE_DIR}/src/MCTruthTrackingAction.cc )

set( headers 
               ${PROJECT_SOURCE_DIR}/include/DetectorConstruction.hh         
               ${PROJECT_SOURCE_DIR}/include/MCTruthEventAction.hh           
               ${PROJECT_SOURCE_DIR}/include/MCTruthTrackInformation.hh      
               ${PROJECT_SOURCE_DIR}/include/PrimaryGeneratorAction.hh
               ${PROJECT_SOURCE_DIR}/include/MCTruthConfig.hh                
               ${PROJECT_SOURCE_DIR}/include/MCTruthManager.hh               
               ${PROJECT_SOURCE_DIR}/include/MCTruthTrackingAction.hh )

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

