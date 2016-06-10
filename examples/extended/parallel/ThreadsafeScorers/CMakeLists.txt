cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

set(name ts_scorers)
project(ThreadsafeScorers)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# macros
file(GLOB macros ${PROJECT_SOURCE_DIR}/*.mac ${PROJECT_SOURCE_DIR}/*.in)

add_executable(${name} ${name}.cc ${headers}
                                  ${sources}
                                  ${macros})

target_link_libraries(${name} ${Geant4_LIBRARIES})

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

