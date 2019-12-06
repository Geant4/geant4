cmake_minimum_required(VERSION 3.1.3 FATAL_ERROR)

set(name ts_scorers)
project(ThreadsafeScorers C CXX)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

option(USE_TIMEMORY "Enable timemory if not built at Geant4 level" OFF)
set(timemory_COMPONENTS "headers;caliper;papi;gotcha;gperftools-cpu;vector"
    CACHE STRING "timemory components")
set(timemory_FIND_COMPONENTS_INTERFACE ${PROJECT_NAME}-timemory)

if(USE_TIMEMORY AND NOT TARGET geant4-timemory)
    find_package(timemory 3.0 REQUIRED COMPONENTS ${timemory_COMPONENTS})
    set(timemory_LIBRARIES ${PROJECT_NAME}-timemory)
endif()

include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# macros
file(GLOB macros ${PROJECT_SOURCE_DIR}/*.mac ${PROJECT_SOURCE_DIR}/*.in)

add_executable(${name} ${name}.cc ${headers}
                                  ${sources}
                                  ${macros})

target_link_libraries(${name} ${Geant4_LIBRARIES} ${timemory_LIBRARIES})

# For IDEs - specifically Xcode
source_group("macros" FILES ${macros})

#==============================================================================#
# - Copy the macros to the build directory
#------------------------------------------------------------------------------#
foreach(_macro ${macros})
    get_filename_component(_macro ${_macro} NAME)
    configure_file(${PROJECT_SOURCE_DIR}/${_macro}
                  ${PROJECT_BINARY_DIR}/${_macro}
                  COPYONLY)
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ${name} DESTINATION bin)
