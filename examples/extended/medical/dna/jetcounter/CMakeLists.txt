cmake_minimum_required(VERSION 3.16...3.21)
project(jetcounter)

# GEANT4 core
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

# root
find_package(ROOT)
include_directories(${ROOT_INCLUDE_DIRS})


set(CMAKE_BUILD_TYPE release)

# User code
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(jetcounter jetcounter.cc ${sources} ${headers})
target_link_libraries(jetcounter ${Geant4_LIBRARIES})
target_link_libraries(jetcounter ${ROOT_LIBRARIES})

set(OTHERS
        init_phys.in primaries.in run.mac vis.mac
        energy.spectrum
        plot.py)

foreach(_other ${OTHERS})
    configure_file(
            ${PROJECT_SOURCE_DIR}/${_other}
            ${PROJECT_BINARY_DIR}/${_other}
            COPYONLY
    )
endforeach()
