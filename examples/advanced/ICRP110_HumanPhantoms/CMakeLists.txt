#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.27)

project(ICRPphantoms)

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
                    ${Geant4_INCLUDE_DIR}
                    ${GDML_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(ICRP110phantoms ICRP110phantoms.cc ${sources} ${headers})
add_executable(ICRP110standalone ICRP110standalone.cc ${sources} ${headers})

target_link_libraries(ICRP110phantoms ${Geant4_LIBRARIES} )
target_link_libraries(ICRP110standalone ${Geant4_LIBRARIES} )

# Depend on data for runtime
add_dependencies(ICRP110phantoms ICRPdata)
add_dependencies(ICRP110standalone ICRPdata)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build human_phantom. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(ICRPphantoms_SCRIPTS
     female.in female_head.in female_trunk.in male.in male_head.in male_trunk.in vis.mac primary.mac ColourMap.dat openGLVis.mac g4views/g4_00.g4view g4views/g4_01.g4view standalone.mac
#    ICRPdata/Data.dat 
#    "ICRPdata/ICRP110_g4dat/AF/*.g4dat"
#    "ICRPdata/ICRP110_g4dat/AM/*.g4dat"
  )
 
foreach(_script ${ICRPphantoms_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#set(ICRPphantoms_DATADIR
#   ICRPdata
#)
        
#file(COPY ${PROJECT_SOURCE_DIR}/${ICRPphantoms_DATADIR} DESTINATION
#     ${PROJECT_BINARY_DIR})

include(ExternalProject)
ExternalProject_Add(ICRPdata
  SOURCE_DIR ${PROJECT_BINARY_DIR}/ICRPdata
  URL https://cern.ch/geant4-data/datasets/examples/advanced/ICRP110Phantoms/ICRPdata.tar.gz
  URL_MD5 2aba5409d4e7b7a2f0328e7d3402f187
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  )

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(ICRPphantoms DEPENDS phantom)
add_custom_target(ICRPstandalone DEPENDS phantom)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ICRP110phantoms DESTINATION bin)
install(TARGETS ICRP110standalone DESTINATION bin)

