#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
project(HepMCEx02)

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
find_package(HepMC3 REQUIRED)
add_definitions(-DG4LIB_USE_HEPMC3)
set (EVENT_RECORD_INCLUDE_DIR ${HEPMC3_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/external)
set (EVENT_RECORD_LIBRARIES ${HEPMC3_LIB})

set (interface_sources
              ${PROJECT_SOURCE_DIR}/src/HepMC3G4AsciiReader.cc
              ${PROJECT_SOURCE_DIR}/src/HepMC3G4AsciiReaderMessenger.cc
              ${PROJECT_SOURCE_DIR}/src/HepMC3G4Interface.cc
              ${PROJECT_SOURCE_DIR}/src/HepMC3G4PythiaInterface.cc
              ${PROJECT_SOURCE_DIR}/src/HepMC3G4PythiaMessenger.cc 
              #${PROJECT_SOURCE_DIR}/src/HEPEvtcom.cc
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
  # ${PROJECT_SOURCE_DIR}/external/upevnt.F
  # ${PROJECT_SOURCE_DIR}/external/upinit.F
  # ${PROJECT_SOURCE_DIR}/external/upveto.F
  ""
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

set (sources 

              ${PROJECT_SOURCE_DIR}/src/H02DetectorConstruction.cc
              ${PROJECT_SOURCE_DIR}/src/H02EventAction.cc
              ${PROJECT_SOURCE_DIR}/src/H02Field.cc
              ${PROJECT_SOURCE_DIR}/src/H02MuonHit.cc
              ${PROJECT_SOURCE_DIR}/src/H02MuonSD.cc
              ${PROJECT_SOURCE_DIR}/src/H02PrimaryGeneratorAction.cc
              ${PROJECT_SOURCE_DIR}/src/H02PrimaryGeneratorMessenger.cc
              ${PROJECT_SOURCE_DIR}/src/H02SteppingAction.cc
) 


set( headers 
                   ${PROJECT_SOURCE_DIR}/include/H02DetectorConstruction.hh      
                   ${PROJECT_SOURCE_DIR}/include/H02MuonSD.hh                    
                   ${PROJECT_SOURCE_DIR}/include/H02EventAction.hh               
                   ${PROJECT_SOURCE_DIR}/include/H02PrimaryGeneratorAction.hh    
                   ${PROJECT_SOURCE_DIR}/include/H02Field.hh                     
                   ${PROJECT_SOURCE_DIR}/include/H02PrimaryGeneratorMessenger.hh 
                   ${PROJECT_SOURCE_DIR}/include/H02MuonHit.hh                   
                   ${PROJECT_SOURCE_DIR}/include/H02SteppingAction.hh            

                   ${PROJECT_SOURCE_DIR}/include/HepMC3G4Interface.hh            
                   ${PROJECT_SOURCE_DIR}/include/HepMCG4AsciiReaderMessenger.hh
                   ${PROJECT_SOURCE_DIR}/include/HepMC3G4AsciiReaderMessenger.hh 
                   ${PROJECT_SOURCE_DIR}/include/HepMCG4AsciiReader.hh           
                   ${PROJECT_SOURCE_DIR}/include/HepMCG4PythiaMessenger.hh                   
                   ${PROJECT_SOURCE_DIR}/include/HepMC3G4PythiaInterface.hh      
                   ${PROJECT_SOURCE_DIR}/include/HepMCG4Interface.hh 
                   ${PROJECT_SOURCE_DIR}/include/HepMC3G4AsciiReader.hh          
                   ${PROJECT_SOURCE_DIR}/include/HepMC3G4PythiaMessenger.hh      
                   ${PROJECT_SOURCE_DIR}/include/HepMCG4PythiaInterface.hh
                   )

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(HepMCEx02 HepMCEx02.cc ${sources} ${generator_sources} ${interface_sources} ${headers})
target_link_libraries(HepMCEx02 ${Geant4_LIBRARIES}
                                ${EVENT_RECORD_LIBRARIES}
                                ${GENERATOR_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build HepMCEx02. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(HepMCEx02_SCRIPTS
    hepmc_ascii.in hepmc_ascii.out hepmc_pygen.in hepmc_pygen.out pyh4mu.mac pyset.mac vis.mac
    data/example_MyPythia.dat
  )

foreach(_script ${HepMCEx02_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS HepMCEx02 DESTINATION bin)

