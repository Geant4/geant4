cmake_minimum_required( VERSION 2.6 FATAL_ERROR )
project( test22_PiPinclusive )
find_package( Geant4 REQUIRED )

include( ${Geant4_USE_FILE} )
include_directories( ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/../src ${PROJECT_SOURCE_DIR}/../include)
file( GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/../src/*.cc)
file( GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/../include/*.hh)

#-----
add_executable( test22_PimP100 PimP100.cc ${sources} ${headers} )
target_link_libraries( test22_PimP100 ${Geant4_LIBRARIES} )

install( TARGETS test22_PimP100 DESTINATION bin )

#-----
add_executable( test22_PimP16 PimP16.cc ${sources} ${headers} )
target_link_libraries( test22_PimP16 ${Geant4_LIBRARIES} )

install( TARGETS test22_PimP16 DESTINATION bin )

#-----
add_executable( test22_PimP40 PimP40.cc ${sources} ${headers} )
target_link_libraries( test22_PimP40 ${Geant4_LIBRARIES} )

install( TARGETS test22_PimP40 DESTINATION bin )

#-----
add_executable( test22_PipP100 PipP100.cc ${sources} ${headers} )
target_link_libraries( test22_PipP100 ${Geant4_LIBRARIES} )

install( TARGETS test22_PipP100 DESTINATION bin )

#-----
add_executable( test22_PipP16 PipP16.cc ${sources} ${headers} )
target_link_libraries( test22_PipP16 ${Geant4_LIBRARIES} )

install( TARGETS test22_PipP16 DESTINATION bin )

#-----
add_executable( test22_PipP175 PipP175.cc ${sources} ${headers} )
target_link_libraries( test22_PipP175 ${Geant4_LIBRARIES} )

install( TARGETS test22_PipP175 DESTINATION bin )

#-----
add_executable( test22_PipP8 PipP8.cc ${sources} ${headers} )
target_link_libraries( test22_PipP8 ${Geant4_LIBRARIES} )

install( TARGETS test22_PipP8 DESTINATION bin )
