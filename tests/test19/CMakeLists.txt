#----------------------------------------------------------------------------
#
#  Note-1: Env.variables G4INSTALL, G4LIB, ROOTSYS need to be setup
#
#  Note-2: Maybe ROOT can be pulled in via ROOT_INCLUDE_DIR - need to double check...
#
# Usage (for local testing/development):
#
# mkdir /tests/test19-build
# cd /tests/test19-build
# cmake -DGeant4_DIR=$G4LIB/Geant4-9.6.0 \
#       ../test19
# make
#
# Note-3: in some cases, one would also use -DGeant4_SOURCE_DIR=$G4INSTALL option,
#         but here's it's not necessary - the system will even tell about that
#


# Specify mandatory
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

#----Included in Geant4------------------------------------------------------------------
#----Standalone--------------------------------------------------------------------------

#------------------------------------------------------------------------------
# - Enforce an out-of-source builds before anything else
#
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(STATUS "Geant4 requires an out-of-source build.")
    message(STATUS "Please remove these files from ${CMAKE_BINARY_DIR} first:")
    message(STATUS "CMakeCache.txt")
    message(STATUS "CMakeFiles")
    message(STATUS "Once these files are removed, create a separate directory")
    message(STATUS "and run CMake from there")
    message(FATAL_ERROR "in-source build detected")
endif()

project(test19)
# message("Project Name: ${CMAKE_PROJECT_NAME}")

#----------------------------------------------------------------------------
# Find Geant4 package, no UI and Vis drivers activated
#
find_package(Geant4 REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# setup include directories
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${PROJECT_SOURCE_DIR}/../test23/shared/g4app/include
                    ${Geant4_INCLUDE_DIR})

#----------------------------------------------------------------------------
#
# sources and headers for this project
#
file(GLOB prj_sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB prj_headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#
#  sources from shared code
#
file(GLOB g4app_sources 
  ${PROJECT_SOURCE_DIR}/../test23/shared/g4app/src/*.cc)
file(GLOB g4app_headers 
  ${PROJECT_SOURCE_DIR}/../test23/shared/g4app/include/*.hh)

#
# find additional ROOT package
#
find_package(ROOT QUIET)
if(NOT ROOT_FOUND)
  message(STATUS "G4 TESTS: ROOT package not found. ")  
  return()
endif()

#
# additional headers and spources from ROOT (if any)
#
include_directories(${PROJECT_SOURCE_DIR}/../test23/shared/rootanalysis/include
		    ${ROOTSYS}/include)
file(GLOB rootana_sources 
  ${PROJECT_SOURCE_DIR}/../test23/shared/rootanalysis/src/*.cc)
file(GLOB rootana_headers 
  ${PROJECT_SOURCE_DIR}/../test23/shared/rootana/include/*.hh)

#
# finally, put all headers and sources together 
#
list(APPEND sources ${prj_sources} ${g4app_sources} ${rootana_sources})
list(APPEND headers ${prj_headers} ${g4app_headers} ${rootana_headers})


#----------------------------------------------------------------------------
#
# specify what executable, and what are the sources for it...
#
add_executable(test19 test19.cc ${sources} ${headers})

#
# ... plus, the external libraries
#
target_link_libraries(test19 ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

# copy over the ctest config
#
set(test19_SCRIPTS ctest19-na61.config)
foreach(_script ${test19_SCRIPTS})
   configure_file(${PROJECT_SOURCE_DIR}/${_script} ${PROJECT_BINARY_DIR}/${_script} COPYONLY)
endforeach()

# alternative approach:
#
# file(GLOB test19_FILES RELATIVE ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/ctest19-na61.config)
# foreach(_file ${test19_FILES})
#    configure_file(
#       ${PROJECT_SOURCE_DIR}/${_file} ${PROJECT_BINARY_DIR}/${_file} COPYONLY
#       )
# endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS test19 DESTINATION bin)

# -------------------------------
# This is a leftover from earlier days, and it's kep here just as an example
#
### if(CMAKE_PROJECT_NAME STREQUAL Geant4)
###  find_package(Geant4 REQUIRED)
###  include(${Geant4_USE_FILE})
###  find_package(ROOT QUIET)
###  if(NOT ROOT_FOUND)
###    message(STATUS "G4 TESTS: ROOT package not found. --> test19 disabled")
###    return()
###  endif()
###  GEANT4_ADD_TEST(test19-build BUILD test19)
####  GEANT4_ADD_TEST(test19-na61 COMMAND test19 ${PROJECT_SOURCE_DIR}/ctest19-na61.config
####                         DEPENDS test19-build
####                         ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
###  GEANT4_ADD_TEST(test19-na61 COMMAND test19 ${CMAKE_CURRENT_SOURCE_DIR}/ctest19-na61.config
###                         DEPENDS test19-build
###                         ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
### endif()
