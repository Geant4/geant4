#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.27)
project(OpNovice)
#add_definitions(-DGEANT4_USE_GDML)
add_compile_definitions(GEANT4_USE_GDML)
#make using gdml the default
OPTION( GEANT4_USE_GDML "Use gdml Option" ON)
# printout the value of GEANT4_USE_GDML
MESSAGE("?GEANT4_USE_GDML " ${GEANT4_USE_GDML})
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
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
if (GEANT4_USE_GDML)
file(GLOB detectors ${PROJECT_SOURCE_DIR}/gdml/*.gdml)
endif()
file(GLOB macros ${PROJECT_SOURCE_DIR}/*.mac)
file(GLOB inputs ${PROJECT_SOURCE_DIR}/*.in)
if (NOT GEANT4_USE_GDML)
list(REMOVE_ITEM sources ${PROJECT_SOURCE_DIR}/src/OpNoviceGDMLDetectorConstructionMessenger.cc)
list(REMOVE_ITEM sources ${PROJECT_SOURCE_DIR}/src/OpNoviceGDMLDetectorConstruction.cc)
list(REMOVE_ITEM headers ${PROJECT_SOURCE_DIR}/include/OpNoviceGDMLDetectorConstructionMessenger.hh)
list(REMOVE_ITEM headers ${PROJECT_SOURCE_DIR}/include/OpNoviceGDMLDetectorConstruction.hh)
endif()
#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(OpNovice OpNovice.cc ${sources} ${headers})
target_link_libraries(OpNovice ${Geant4_LIBRARIES} )

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build OpNovice. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(OpNovice_SCRIPTS
    OpNovice.out
    OpNovice.in
    optPhoton.mac
    gui.mac
    vis.mac
  )

foreach(_script ${OpNovice_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS OpNovice DESTINATION bin)
if (GEANT4_USE_GDML)
install(FILES ${detectors} ${macros} ${inputs} DESTINATION bin)
else()
install(FILES ${macros} ${inputs} DESTINATION bin)
endif()