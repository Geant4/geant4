# - CmakeLists.txt for G4mpi interface

#------------------------------------------------------------------------------
set(_projname libG4mpi)
set(_targetname G4mpi)
project(${_projname})

#------------------------------------------------------------------------------
#Manadatory dependencies

cmake_minimum_required(VERSION 3.16...3.21)
find_package(MPI REQUIRED)
find_package(Geant4 10.2.0 REQUIRED)
include(${Geant4_USE_FILE})

#------------------------------------------------------------------------------
message("-- G4 Examples: ${_projname} uses includes from: "
        "${MPI_CXX_INCLUDE_PATH}")
message("-- G4 Examples: ${_projname} uses libraries: ${MPI_CXX_LIBRARIES}")

#Version of this package is the same as the G4 version
set(${_targetname}_VERSION ${Geant4_VERSION})
set(${_targetname}_VERSION_MAJOR ${Geant4_VERSION_MAJOR})
set(${_targetname}_VERSION_MINOR ${Geant4_VERSION_MINOR})
set(${_targetname}_VERSION_PATCH ${Geant4_VERSION_PATCH})

#------------------------------------------------------------------------------
#Options for this package
option(BUILD_SHARED_LIBS "If true build shared library" ON)
option(BUILD_STATIC_LIBS "If true build static library" OFF)
option(G4MPI_OLD_MPI "If true use old signatures for MPI_[Un]Pack functions" 
       OFF)
mark_as_advanced(G4MPI_OLD_MPI)

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH
  "Installation directory for header files")
#INSTALL_CMAKE_DIR is set later
 
#------------------------------------------------------------------------------
#Setup compilation options specific to this project
add_definitions(-DTOOLS_USE_NATIVE_MPI) 
if(G4MPI_OLD_MPI)
  message("-- G4 Examples: ${_projname} using old non-const signatures for" 
          "MPI_[Un]Pack functions")
  add_definitions(-DTOOLS_USE_MPI_PACK_NOT_CONST)
  add_definitions(-DG4MPI_USE_MPI_PACK_NOT_CONST)
  # Prevent warning with c++11 on old MPI 
  add_definitions(-Wno-literal-suffix)
endif()

#------------------------------------------------------------------------------
#TODO: ??? What is G4 policy on this?
#if(Geant4_static_FOUND)
#  set(BUILD_STATIC_LIBS ON)
#  set(BUILD_SHARED_LIBS OFF)
#else()
#  set(BUILD_STATIC_LIBS OFF)
#  set(BUILD_SHARED_LIBS ON)
#endif()

#------------------------------------------------------------------------------
# Define library
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include 
                    ${CMAKE_CURRENT_SOURCE_DIR}/analysis/include 
                    ${Geant4_INCLUDE_DIR} 
                    ${MPI_CXX_INCLUDE_PATH})
link_directories(${MPI_CXX_LIBRARY_DIRS})

set(_sources 
  src/G4MPIbatch.cc
  src/G4MPIextraWorker.cc
  src/G4MPImanager.cc
  src/G4MPImessenger.cc
  src/G4MPIrandomSeedGenerator.cc
  src/G4MPIsession.cc
  src/G4MPIstatus.cc
  src/G4UImpish.cc
  src/G4VMPIseedGenerator.cc
  src/G4VMPIsession.cc
  src/G4MPIscorerMerger.cc
  src/G4MPIhistoMerger.cc
  src/G4MPIntupleMerger.cc
  src/G4VUserMPIrunMerger.cc
  src/G4MPIutils.cc
  analysis/src/G4RootMpiAnalysisManager.cc
  analysis/src/G4RootMpiNtupleFileManager.cc
  analysis/src/G4RootMpiNtupleManager.cc
  analysis/src/G4RootMpiPNtupleManager.cc
)

#TODO: This works only if both BUILD_STATIC_LIBS=ON && BUILD_SHARED_LIBS=OFF
#      are explicitly specified, what should be the behavior if only
#      BUILD_STATIC_LIBS is specified?
#      In Geant4 both are created: libXXX.so and libXXX.a with two target names
#      XXX and XXX-static. I should study how to do that...
if(NOT BUILD_SHARED_LIBS AND NOT BUILD_STATIC_LIBS)
    message(FATAL "-- G4 Examples: ${_projname} neither static of shared build selected")
    return()
endif()
if(BUILD_STATIC_LIBS)
  add_library(${_targetname}-static STATIC ${_sources})
  set_target_properties(${_targetname}-static PROPERTIES OUTPUT_NAME ${_targetname})
  target_link_libraries(${_targetname}-static ${MPI_CXX_LIBRARIES} ${Geant4_LIBRARIES})
endif()
if(BUILD_SHARED_LIBS)
  add_library(${_targetname} SHARED ${_sources})
  target_link_libraries(${_targetname} ${MPI_CXX_LIBRARIES} ${Geant4_LIBRARIES})
endif()


# headers
set(HEADERS
  include/G4MPIbatch.hh
  include/G4MPIextraWorker.hh
  include/G4MPImanager.hh
  include/G4MPImessenger.hh
  include/G4MPIrandomSeedGenerator.hh
  include/G4MPIsession.hh
  include/G4MPIstatus.hh
  include/G4UImpish.hh
  include/G4VMPIextraWorker.hh
  include/G4VMPIseedGenerator.hh
  include/G4VMPIsession.hh
  include/G4MPIscorerMerger.hh
  include/G4MPIrunMerger.hh
  include/G4MPIhistoMerger.hh
  include/G4MPIntupleMerger.hh
  include/G4VUserMPIrunMerger.hh
  include/G4MPIutils.hh
  analysis/include/G4RootMpiAnalysisManager.hh
  analysis/include/G4RootMpiNtupleFileManager.hh
  analysis/include/G4RootMpiNtupleManager.hh
  analysis/include/G4RootMpiPNtupleDescription.hh
  analysis/include/G4RootMpiPNtupleManager.hh

)

#------------------------------------------------------------------------------
# Overwrite "lib" path w/ "lib64" if needed
set(_LIBDIR_DEFAULT "lib")
set(_dolib64 FALSE)
if(NOT DEFINED CMAKE_INSTALL_LIBDIR)
  set(_dolib64 TRUE)
endif()
#If this is built as part of G4 CMAKE_INSTALL_LIBDIR is defined, we need to
#so force check of what we should use
include(G4DeveloperAPI OPTIONAL RESULT_VARIABLE _internal_build)
if(_internal_build)
  set(_dolib64 TRUE)
endif()
if(${_dolib64})
  # Override this default 'lib' with 'lib64' iff:
  #  - we are on Linux system but NOT cross-compiling
  #  - we are NOT on debian
  #  - we are on a 64 bits system
  # reason is: amd64 ABI: http://www.x86-64.org/documentation/abi.pdf
  # Note that the future of multi-arch handling may be even
  # more complicated than that: http://wiki.debian.org/Multiarch
  if(CMAKE_SYSTEM_NAME MATCHES "Linux"
      AND NOT CMAKE_CROSSCOMPILING
      AND NOT EXISTS "/etc/debian_version")
    if(NOT DEFINED CMAKE_SIZEOF_VOID_P)
      message(AUTHOR_WARNING
        "Unable to determine default CMAKE_INSTALL_LIBDIR directory because "
        "no target architecture is known. "
        "Please enable at least one language before including GNUInstallDirs.")
    else()
      if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
        set(_LIBDIR_DEFAULT "lib64")
      endif()
    endif()
  endif()
  set(CMAKE_INSTALL_LIBDIR "${_LIBDIR_DEFAULT}" CACHE PATH "object code libraries (${_LIBDIR_DEFAULT})")
  set(INSTALL_LIB_DIR ${CMAKE_INSTALL_LIBDIR})
endif()

#Set Location of .cmake files
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR CMake)
else()
  set(DEF_INSTALL_CMAKE_DIR 
        ${_LIBDIR_DEFAULT}/${_targetname}-${${_targetname}_VERSION})
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
  "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

#------------------------------------------------------------------------------
#install
if(BUILD_SHARED_LIBS)
  INSTALL(TARGETS ${_targetname}
          EXPORT ${_targetname}Targets 
          LIBRARY DESTINATION ${INSTALL_LIB_DIR} COMPONENT shlib)
endif()
if(BUILD_STATIC_LIBS)
   INSTALL(TARGETS ${_targetname}-static
           EXPORT ${_targetname}Targets 
           ARCHIVE DESTINATION ${INSTALL_LIB_DIR} COMPONENT dev)
endif()
INSTALL(FILES ${HEADERS} DESTINATION include)

#G4mpiTargets.cmake files
if(BUILD_SHARED_LIBS)
    export(TARGETS ${_targetname}
           FILE "${PROJECT_BINARY_DIR}/${_targetname}Targets.cmake")
endif()
if(BUILD_STATIC_LIBS)
    export(TARGETS ${_targetname}-static
           FILE "${PROJECT_BINARY_DIR}/${_targetname}Targets.cmake")
endif()
# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE ${_projname})

# Create the G4mpiConfig.cmake and G4mpiConfigVersion files
set(TARGET_STATIC OFF)
if(BUILD_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)#If both are activated, prefer shared
   set(TARGET_STATIC ON)
endif()
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
       "${INSTALL_INCLUDE_DIR}")
# ... for the build tree
set(TARGET_NAME ${_targetname})
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/include")
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
   "${INSTALL_INCLUDE_DIR}")
configure_file(G4mpiConfig.cmake.in
  "${PROJECT_BINARY_DIR}/${_targetname}Config.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${${_targetname}_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(${_targetname}Config.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${_targetname}Config.cmake" @ONLY)
# ... for both
set(THE_VERSION ${${_targetname}_VERSION})
configure_file(G4mpiConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}/${_targetname}ConfigVersion.cmake" @ONLY)
  
# Install the G4mpiConfig.cmake 
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${_targetname}Config.cmake"
  "${PROJECT_BINARY_DIR}/${_targetname}ConfigVersion.cmake" 
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
# Install the export set for use with the install-tree
install(EXPORT ${_targetname}Targets DESTINATION
  "${INSTALL_CMAKE_DIR}" COMPONENT dev)


