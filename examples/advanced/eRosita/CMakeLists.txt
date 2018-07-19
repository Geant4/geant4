#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(eRosita)

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


#-----------------------------------------------------------------------------------
# Notice G4PIIDATA environment variable as appropriate to get acces to the PIXE data 
#
# message("G4PIIDATA: " $ENV{G4PIIDATA})
if(DEFINED ENV{G4PIIDATA})
  message(STATUS "G4PIIDATA environment variable defined")
else()
  message(STATUS "Please set G4PIIDATA enviroment (e.g. G4PII1.3) in your system!")    
endif()


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/application/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/application/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/application/include/*.hh)

#-----------------------------------------------------------------------------------
# Setup Geant4 physics mode 
#
# ON : simulate Geant4 with PIXE data
# OFF: simulate Geant4 with original Geant4 low energy physics processes
#       based on EEDL/EPDL/EADL data libraries
option(PHYSICS_MODE "Please select Geant4 Physics model. ON:PIXE, OFF:OriginalGeant4Physics" ON)
if(PHYSICS_MODE)
  add_executable(eRosita eRosita.cc ${sources} ${headers})
else()
  add_definitions(-DORIGINAL_PHYSICS)
  file(GLOB sourcesOriginal ${PROJECT_SOURCE_DIR}/physics/src/*.cc)
  file(GLOB headersOriginal ${PROJECT_SOURCE_DIR}/physics/include/*.hh)
  include_directories(${PROJECT_SOURCE_DIR}/physics/include)
  add_executable(eRosita eRosita.cc ${sources} ${headers} 
                                  ${sourcesOriginal} ${headersOriginal} )
endif()


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
target_link_libraries(eRosita ${Geant4_LIBRARIES})


#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build eRosita. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEeRosita_SCRIPTS
  eRosita.in
  eRositaLong.in
)

foreach(_script ${EXAMPLEeRosita_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
#add_custom_target(eRosita DEPENDS eRosita)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS eRosita DESTINATION bin)
