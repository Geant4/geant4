cmake_minimum_required( VERSION 2.6 FATAL_ERROR )
project( test22_PipPchan )
find_package( Geant4 REQUIRED )

include( ${Geant4_USE_FILE} )
include_directories( ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/../src ${PROJECT_SOURCE_DIR}/../include)
file( GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/../src/*.cc)
file( GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/../include/*.hh)

#-----
add_executable( test22_PipPchan PipPchan.cc ${sources} ${headers} )
target_link_libraries( test22_PipPchan ${Geant4_LIBRARIES} )

install( TARGETS test22_PipPchan DESTINATION bin )
