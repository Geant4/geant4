#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(HepMCEx01)

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
# Find HepMC3 (required package)
#

add_definitions(-DG4LIB_USE_HEPMC3)
find_package(HepMC3 REQUIRED)
set (EVENT_RECORD_INCLUDE_DIR ${HEPMC3_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/external)
set (EVENT_RECORD_LIBRARIES ${HEPMC3_LIB})
set (interface_sources
                        ${PROJECT_SOURCE_DIR}/src/HepMC3G4AsciiReader.cc
                        ${PROJECT_SOURCE_DIR}/src/HepMC3G4AsciiReaderMessenger.cc
                        ${PROJECT_SOURCE_DIR}/src/HepMC3G4Interface.cc
                        ${PROJECT_SOURCE_DIR}/src/HepMC3G4PythiaInterface.cc
                        ${PROJECT_SOURCE_DIR}/src/HepMC3G4PythiaMessenger.cc 
                        )

#----------------------------------------------------------------------------
# Find Pythia6 (optional package)
#
find_package(Pythia6)
include(CheckLanguage)
check_language(Fortran)
if(PYTHIA6_FOUND AND CMAKE_Fortran_COMPILER)
  message(STATUS "G4 Examples: Pythia6 found. --> HepMCEx01 example with Pythia6 enabled.") 
  add_definitions(-DG4LIB_USE_PYTHIA) 
  enable_language(Fortran)
   message(STATUS "Fortran compiler found  Pythia6  were found ${Fortran_COMPILER_NAME} ")
   set( GENERATOR_LIBRARIES ${PYTHIA6_LIBRARIES} LHAPDF gfortran)
   set (generator_sources 
   #${PROJECT_SOURCE_DIR}/external/upevnt.F
   #${PROJECT_SOURCE_DIR}/external/upinit.F
   #${PROJECT_SOURCE_DIR}/external/upveto.F
   )

else()
  set(GENERATOR_LIBRARIES "")  
 set (generator_sources "")
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
                    ${EVENT_RECORD_INCLUDE_DIR})
set(sources 
              ${PROJECT_SOURCE_DIR}/src/ExN04CalorimeterHit.cc
              ${PROJECT_SOURCE_DIR}/src/ExN04CalorimeterParametrisation.cc
              ${PROJECT_SOURCE_DIR}/src/ExN04CalorimeterROGeometry.cc
              ${PROJECT_SOURCE_DIR}/src/ExN04CalorimeterSD.cc
              ${PROJECT_SOURCE_DIR}/src/ExN04DetectorConstruction.cc
              ${PROJECT_SOURCE_DIR}/src/ExN04DetectorParameterDef.icc
              ${PROJECT_SOURCE_DIR}/src/ExN04EventAction.cc
              ${PROJECT_SOURCE_DIR}/src/ExN04Field.cc
              ${PROJECT_SOURCE_DIR}/src/ExN04MuonHit.cc
              ${PROJECT_SOURCE_DIR}/src/ExN04MuonSD.cc
              ${PROJECT_SOURCE_DIR}/src/ExN04PrimaryGeneratorAction.cc
              ${PROJECT_SOURCE_DIR}/src/ExN04PrimaryGeneratorMessenger.cc
              ${PROJECT_SOURCE_DIR}/src/ExN04RunAction.cc
              ${PROJECT_SOURCE_DIR}/src/ExN04StackingAction.cc
              ${PROJECT_SOURCE_DIR}/src/ExN04StackingActionMessenger.cc
              ${PROJECT_SOURCE_DIR}/src/ExN04SteppingAction.cc
              ${PROJECT_SOURCE_DIR}/src/ExN04SteppingVerbose.cc
              ${PROJECT_SOURCE_DIR}/src/ExN04TrackerHit.cc
              ${PROJECT_SOURCE_DIR}/src/ExN04TrackerParametrisation.cc
              ${PROJECT_SOURCE_DIR}/src/ExN04TrackerSD.cc
              ${PROJECT_SOURCE_DIR}/src/ExN04TrackingAction.cc
)


set( headers 
                    ${PROJECT_SOURCE_DIR}/include/ExN04CalorimeterHit.hh                  
                    ${PROJECT_SOURCE_DIR}/include/ExN04Field.hh                           
                    ${PROJECT_SOURCE_DIR}/include/ExN04SteppingAction.hh                  
                    ${PROJECT_SOURCE_DIR}/include/ExN04CalorimeterParametrisation.hh      
                    ${PROJECT_SOURCE_DIR}/include/ExN04MuonHit.hh                         
                    ${PROJECT_SOURCE_DIR}/include/ExN04SteppingVerbose.hh                 
                    ${PROJECT_SOURCE_DIR}/include/ExN04CalorimeterROGeometry.hh           
                    ${PROJECT_SOURCE_DIR}/include/ExN04MuonSD.hh                          
                    ${PROJECT_SOURCE_DIR}/include/ExN04TrackerHit.hh                      
                    ${PROJECT_SOURCE_DIR}/include/ExN04CalorimeterSD.hh                   
                    ${PROJECT_SOURCE_DIR}/include/ExN04PrimaryGeneratorAction.hh          
                    ${PROJECT_SOURCE_DIR}/include/ExN04TrackerParametrisation.hh          
                    ${PROJECT_SOURCE_DIR}/include/ExN04DetectorConstruction.hh            
                    ${PROJECT_SOURCE_DIR}/include/ExN04PrimaryGeneratorMessenger.hh       
                    ${PROJECT_SOURCE_DIR}/include/ExN04TrackerSD.hh                       
                    ${PROJECT_SOURCE_DIR}/include/ExN04DetectorParameterDef.hh            
                    ${PROJECT_SOURCE_DIR}/include/ExN04RunAction.hh                       
                    ${PROJECT_SOURCE_DIR}/include/ExN04TrackingAction.hh                  
                    ${PROJECT_SOURCE_DIR}/include/ExN04DummySD.hh                         
                    ${PROJECT_SOURCE_DIR}/include/ExN04StackingAction.hh                  
                    ${PROJECT_SOURCE_DIR}/include/ExN04EventAction.hh                     
                    ${PROJECT_SOURCE_DIR}/include/ExN04StackingActionMessenger.hh         

                    ${PROJECT_SOURCE_DIR}/include/HepMCG4Interface.hh
                    ${PROJECT_SOURCE_DIR}/include/HepMCG4AsciiReaderMessenger.hh
                    ${PROJECT_SOURCE_DIR}/include/HepMC3G4AsciiReader.hh                  
                    ${PROJECT_SOURCE_DIR}/include/HepMCG4PythiaInterface.hh
                    ${PROJECT_SOURCE_DIR}/include/HepMC3G4AsciiReaderMessenger.hh         
                    ${PROJECT_SOURCE_DIR}/include/HepMCG4PythiaMessenger.hh
                    ${PROJECT_SOURCE_DIR}/include/HepMCG4AsciiReader.hh
                    ${PROJECT_SOURCE_DIR}/include/HepMC3G4PythiaMessenger.hh
                    ${PROJECT_SOURCE_DIR}/include/HepMC3G4PythiaInterface.hh
                    ${PROJECT_SOURCE_DIR}/include/HepMC3G4Interface.hh 
                    )

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(HepMCEx01 HepMCEx01.cc ${sources} ${generator_sources} ${interface_sources} ${headers})
target_link_libraries(HepMCEx01 ${Geant4_LIBRARIES}
                                ${EVENT_RECORD_LIBRARIES}
                                ${GENERATOR_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build HepMCEx01. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(HepMCEx01_SCRIPTS
    hepmc_ascii.in hepmc_ascii.out hepmc_pygen.in hepmc_pygen.out pyh4mu.mac pyset.mac vis.mac
    data/example_MyPythia.dat
  )

foreach(_script ${HepMCEx01_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS HepMCEx01 DESTINATION bin)

