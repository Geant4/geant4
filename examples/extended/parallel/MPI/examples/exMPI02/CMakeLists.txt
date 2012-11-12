# - CmakeLists.txt for building an application

#------------------------------------------------------------------------------
project(exMPI02)
#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

SET(CMAKE_MODULE_PATH
    ${PROJECT_SOURCE_DIR}
    ${CMAKE_MODULE_PATH})

# check MPI package...
find_package(MPI REQUIRED)
# modify these var. if needed
set(CMAKE_CXX_COMPILER mpicxx)
#set(CMAKE_CXX_INCLUDE_PATH )

# check ROOT package...
find_package(ROOT REQUIRED)

#------------------------------------------------------------------------------
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

#------------------------------------------------------------------------------
include(${Geant4_USE_FILE})

SET(Geant4_LIBRARIES ${Geant4_LIBRARIES} G4mpi)
ADD_LIBRARY(G4mpi SHARED IMPORTED)
SET_PROPERTY(TARGET G4mpi APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
SET_TARGET_PROPERTIES(G4mpi PROPERTIES
  IMPORTED_LOCATION_RELEASE "${Geant4_DIR}/../libG4mpi.so"
  IMPORTED_SONAME_RELEASE "libG4mpi.so"
)

#include_directories(${PROJECT_SOURCE_DIR}/include)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIR})


add_executable(exMPI02 exMPI02.cc ${sources} ${headers})
target_link_libraries(exMPI02 ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
set_target_properties(exMPI02
                      PROPERTIES INSTALL_RPATH
                      ${Geant4_DIR}/..
                      BUILD_WITH_INSTALL_RPATH TRUE)

install(TARGETS exMPI02 RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
