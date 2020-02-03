#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(medical_linac)

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

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(ml2 ml2.cc ${sources} ${headers})
target_link_libraries(ml2 ${Geant4_LIBRARIES} )

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build medical_linac. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(medical_linac_SCRIPTS
    ml2.mac  
    batch.mac
    acc1.mac
    macroAndData/macro_files/accelerator/acc1.mac
    macroAndData/macro_files/accelerator/acc2_6MeV_3x3.mac
    macroAndData/macro_files/accelerator/acc2_6MeV_10x10.mac
    macroAndData/macro_files/accelerator/acc2_15MeV_3x3.mac
    macroAndData/macro_files/accelerator/accSaturn.mac
    macroAndData/macro_files/main/batch.mac 
    macroAndData/macro_files/main/ml2.mac 
    macroAndData/macro_files/main/vis.mac
    macroAndData/experimental/VoxelsTest_In.txt
    macroAndData/experimental/EXP_Saturn_Cal_Data.txt
    macroAndData/results/PhSp/PhSp_Acc1.txt
    macroAndData/results/ROG/ROG_Out.txt
    macroAndData/results/experimental/VoxelsTest_Out.m   
  )

foreach(_script ${medical_linac_SCRIPTS})
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
add_custom_target(medical_linac DEPENDS ml2)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ml2 DESTINATION bin)

