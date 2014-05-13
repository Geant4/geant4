#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(go)

set(EXT_INC_DIR)
set(EXT_LIBS)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available Vis drivers by default
# You can set WITH_GEANT4_VIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_VIS "Build example with Geant4 Vis drivers" ON)
if(WITH_GEANT4_VIS)
  find_package(Geant4 REQUIRED vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

list(APPEND EXT_INC_DIR ${Geant4_INCLUDE_DIR})
list(APPEND EXT_LIBS ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Setup ROOT
# Au cas ou tu veuilles compiler avec root plus tard ...
option(useROOT "Build with ROOT" OFF)

message(STATUS "useROOT is ${useROOT}")

if(useROOT)
    find_package(ROOT REQUIRED)
    if(ROOT_FOUND)
    include_directories(${ROOT_INCLUDE_DIR})
    message(STATUS ROOT_INCLUDE_DIR = ${ROOT_INCLUDE_DIR})
    message(STATUS ROOT_LIBRARIES = ${ROOT_LIBRARIES})
    set(CMAKE_CXX_FLAGS "-DUSE_ROOT")
    add_definitions(-DUSE_ROOT)
    LIST(APPEND EXTRA_ROOT_LIBRARIES ${ROOT_LIBRARIES}) #-lTreePlayer -lProof -lProofDraw -lProofPlayer  -lTreeViewer -lGui)
    list(APPEND EXT_INC_DIR ${ROOT_INCLUDE_DIR})
    list(APPEND EXT_LIBS ${EXTRA_ROOT_LIBRARIES})    
    endif(ROOT_FOUND)
endif(useROOT)

file(GLOB tests_files ${PROJECT_SOURCE_DIR}/*.cc)

foreach(_tests ${tests_files})
 include(${Geant4_USE_FILE})
 include_directories(${EXT_INC_DIR})
 include_directories(${PROJECT_SOURCE_DIR}/include
		    ${EXT_INC_DIR})
 GET_FILENAME_COMPONENT(tests_files_we ${_tests} NAME_WE)
 add_executable(${tests_files_we} ${_tests})
 target_link_libraries(${tests_files_we} ${Geant4_LIBRARIES})
endforeach()
