
#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.16...3.21)
project(plotG)

#----------------------------------------------------------------------------
# ROOT
#
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS Gui)
add_definitions(-DUSE_ROOT)

# Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

#----------------------------------------------------------------------------
# Setup include directory for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#

EXECUTE_PROCESS(COMMAND ${ROOT_rootcint_CMD} -f ${PROJECT_BINARY_DIR}/CanvasInTab_Dict.cc
                                                ${PROJECT_SOURCE_DIR}/include/CanvasInTab.hh
                                                ${PROJECT_SOURCE_DIR}/include/CanvasInTab_LinkDef.hh)

file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc
                  ${PROJECT_BINARY_DIR}/CanvasInTab_Dict.cc)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#

add_executable(plotG plotG.cc ${sources} ${headers})
target_link_libraries(plotG ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#

install(TARGETS plotG DESTINATION bin)
