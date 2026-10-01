#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.27)
project(Hadr02)

#----------------------------------------------------------------------------
# Prepend our own CMake Modules to the search path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})

#----------------------------------------------------------------------------
# Check whether CRMC should be used or not

option(G4_USE_CRMC "Using CRMC" OFF) 
if(G4_USE_CRMC)
  message(STATUS "G4_USE_CRMC ON : Using CRMC for building Hadr02")
  add_definitions(-DG4_USE_CRMC)
  message(STATUS "CMAKE_MODULE_PATH =" ${CMAKE_MODULE_PATH} " ; PROJECT_SOURCE_DIR = " ${PROJECT_SOURCE_DIR})
  find_package(CRMC REQUIRED)
  if(CRMC_FOUND)
    message(STATUS "CRMC was found :" ${CRMC_INCLUDE_DIR})
  else()
    message(FATAL_ERROR "CRMC was NOT found !")
  endif()
else()
  message(STATUS "G4_USE_CRMC OFF : NOT using CRMC for building Hadr02")
endif()

#----------------------------------------------------------------------------
# UrQMD 1.3cr interface (OFF by default). Requires the external, licensed
# UrQMD-1.3cr FORTRAN sources unpacked in urqmd1_3/urqmd-1.3cr/. CMake compiles
# them (plus the provided urqmd1_3/g4urqmdblockdata.f) into a static library;
# the Fortran runtime is resolved via CMake's Fortran support.
option(G4_USE_URQMD "Build the UrQMD 1.3cr interface" OFF)
if(G4_USE_URQMD)
  message(STATUS "G4_USE_URQMD ON : building UrQMD 1.3cr for Hadr02")
  set(URQMD_DIR "${PROJECT_SOURCE_DIR}/urqmd1_3/urqmd-1.3cr")
  if(NOT EXISTS "${URQMD_DIR}/UrqmdLink.f")
    message(FATAL_ERROR
      "UrQMD 1.3cr sources not available in ${URQMD_DIR} - "
      "see README.md for details on obtaining them")
  endif()
  enable_language(Fortran)
  add_definitions(-DG4_USE_URQMD)
  # UrQMD physics sources (GNUmakefile SRC set, minus UrqmdMain.f, whose main()
  # clashes with Geant4's) from the unpacked tarball ...
  set(URQMD_FSRC
    UrqmdLink.f addpart.f angdis.f anndec.f blockres.f boxprg.f cascinit.f coload.f
    dectim.f delpart.f detbal.f dwidth.f error.f getmass.f getspin.f gnuranf.f init.f
    input.f iso.f ityp2pdg.f jdecay2.f make22.f numrec.f output.f paulibl.f proppot.f
    saveinfo.f scatter.f siglookup.f string.f tabinit.f urqmd.f whichres.f)
  list(TRANSFORM URQMD_FSRC PREPEND ${URQMD_DIR}/)
  # ... plus the Geant4-provided block data that ships with the example.
  list(APPEND URQMD_FSRC ${PROJECT_SOURCE_DIR}/urqmd1_3/g4urqmdblockdata.f)
  add_library(urqmd1_3 STATIC ${URQMD_FSRC})
  set_target_properties(urqmd1_3 PROPERTIES POSITION_INDEPENDENT_CODE ON)
  # legacy fixed-form FORTRAN; the type/name-mangling flags must match the
  # extern "C" declarations in G4UrQMD1_3Interface.hh (gfortran assumed).
  target_compile_options(urqmd1_3 PRIVATE
    -fno-automatic -fdefault-real-8 -fdefault-double-8 -ff2c -fno-second-underscore
    -ffixed-line-length-132 -frepack-arrays -std=legacy -fallow-argument-mismatch)
endif()

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
include_directories(${CRMC_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${CRMC_INCLUDE_DIR}/*.h ${PROJECT_SOURCE_DIR}/include/*.hh)

# Compile the guarded UrQMD C++ interface even when the licensed FORTRAN
# sources are not available. Its objects are not added to Hadr02; the target
# only keeps the interface checked against changes in the Geant4 C++ API.
if(NOT G4_USE_URQMD)
  set(URQMD_INTERFACE_SOURCES
    ${PROJECT_SOURCE_DIR}/src/G4UrQMD1_3Model.cc
    ${PROJECT_SOURCE_DIR}/src/HadronPhysicsUrQMD.cc
    ${PROJECT_SOURCE_DIR}/src/HistoManager.cc
    ${PROJECT_SOURCE_DIR}/src/IonUrQMDPhysics.cc
    ${PROJECT_SOURCE_DIR}/src/UrQMD.cc
    ${PROJECT_SOURCE_DIR}/src/UrQMDAntiBarionBuilder.cc
    ${PROJECT_SOURCE_DIR}/src/UrQMDNeutronBuilder.cc
    ${PROJECT_SOURCE_DIR}/src/UrQMDPiKBuilder.cc
    ${PROJECT_SOURCE_DIR}/src/UrQMDProtonBuilder.cc)
  add_library(Hadr02UrQMDCompileCheck OBJECT ${URQMD_INTERFACE_SOURCES})
  target_compile_definitions(Hadr02UrQMDCompileCheck PRIVATE G4_USE_URQMD)
  target_include_directories(Hadr02UrQMDCompileCheck PRIVATE ${PROJECT_SOURCE_DIR}/include)
  # Propagate Geant4's compile requirements. OBJECT libraries have no link step.
  target_link_libraries(Hadr02UrQMDCompileCheck PRIVATE ${Geant4_LIBRARIES})
endif()

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(Hadr02 Hadr02.cc ${sources} ${headers})
target_link_libraries(Hadr02 ${CRMC_LIBRARIES} ${Geant4_LIBRARIES})
if(TARGET Hadr02UrQMDCompileCheck)
  add_dependencies(Hadr02 Hadr02UrQMDCompileCheck)
endif()
if(G4_USE_URQMD)
  # CMake propagates the Fortran runtime automatically when linking a Fortran
  # static library into a C++ target, so no explicit gfortran path is needed.
  target_link_libraries(Hadr02 urqmd1_3)
endif()

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build Hadr02. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(Hadr02_SCRIPTS
    hadr02.in hadr02.out urqmd.in crmc.in vis.mac
  )

foreach(_script ${Hadr02_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

if(G4_USE_CRMC)
  if(CRMC_FOUND)
    configure_file(${CRMC_BUILD_DIR}/crmc.param ${PROJECT_BINARY_DIR} COPYONLY)
  endif()
endif()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS Hadr02 DESTINATION bin)

