cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

set(name exampleP01)
project(${name})
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

find_package(ROOT QUIET)
if(NOT ROOT_FOUND)
  message(STATUS "G4 Examples: ROOT package not found. --> Example P01 disabled")
  return()
endif()


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIR})
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)

REFLEX_GENERATE_DICTIONARY(Classes include/Classes.hh SELECTION xml/selection.xml)
add_library(ClassesDict SHARED ${sources} Classes.cpp)

set(libsuffix .so)

set(ROOT_LIBRARY_PROPERTIES ${ROOT_LIBRARY_PROPERTIES}
      SUFFIX ${libsuffix})

set_target_properties(ClassesDict PROPERTIES ${ROOT_LIBRARY_PROPERTIES})

target_link_libraries(ClassesDict  ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} Cintex Reflex)

add_executable(${name} ${name}.cc )
target_link_libraries(${name} ClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} Cintex)

add_executable(readHits readHits.cc )
target_link_libraries(readHits ClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} Cintex)

