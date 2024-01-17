#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.16...3.27)
project(runAna)
set(CMAKE_CXX_STANDARD 17)

# Find ROOT package
find_package(ROOT REQUIRED)
#----------------------------------------------------------------------------
# Setup ROOT include directories and compile definitions
# Setup include directory for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${ROOT_INCLUDE_DIRS}  
                    ${PROJECT_SOURCE_DIR}/dnadamage/include  
                    ${PROJECT_SOURCE_DIR}/repairmodels/include        
                    )

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc
    ${PROJECT_SOURCE_DIR}/dnadamage/src/*.cc
    ${PROJECT_SOURCE_DIR}/repairmodels/src/*.cc  
    )
		   
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh
    ${PROJECT_SOURCE_DIR}/dnadamage/include/*.hh
    ${PROJECT_SOURCE_DIR}/repairs/include/*.hh
    )

add_executable(${PROJECT_NAME} main.cc ${sources} ${headers})
target_link_libraries(${PROJECT_NAME} ${ROOT_LIBRARIES} )
