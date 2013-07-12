cmake_minimum_required( VERSION 2.6 FATAL_ERROR )
project( test22_PbarPinclusive )
find_package( Geant4 REQUIRED )

include( ${Geant4_USE_FILE} )
include_directories( ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/../src ${PROJECT_SOURCE_DIR}/../include)
file( GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/../src/*.cc)
file( GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/../include/*.hh)

#-----
add_executable( test22_PbarP100 PbarP100.cc ${sources} ${headers} )
target_link_libraries( test22_PbarP100 ${Geant4_LIBRARIES} )

install( TARGETS test22_PbarP100 DESTINATION bin )

#-----
add_executable( test22_PbarP14p8 PbarP14p8.cc ${sources} ${headers} )
target_link_libraries( test22_PbarP14p8 ${Geant4_LIBRARIES} )

install( TARGETS test22_PbarP14p8 DESTINATION bin )

#-----
add_executable( test22_PbarP22p4 PbarP22p4.cc ${sources} ${headers} )
target_link_libraries( test22_PbarP22p4 ${Geant4_LIBRARIES} )

install( TARGETS test22_PbarP22p4 DESTINATION bin )

#-----
add_executable( test22_PbarP32 PbarP32.cc ${sources} ${headers} )
target_link_libraries( test22_PbarP32 ${Geant4_LIBRARIES} )

install( TARGETS test22_PbarP32 DESTINATION bin )

#-----
add_executable( test22_PbarP4p6 PbarP4p6.cc ${sources} ${headers} )
target_link_libraries( test22_PbarP4p6 ${Geant4_LIBRARIES} )

install( TARGETS test22_PbarP4p6 DESTINATION bin )

#-----
add_executable( test22_PbarP5p7 PbarP5p7.cc ${sources} ${headers} )
target_link_libraries( test22_PbarP5p7 ${Geant4_LIBRARIES} )

install( TARGETS test22_PbarP5p7 DESTINATION bin )

#-----
add_executable( test22_PbarP7p3 PbarP7p3.cc ${sources} ${headers} )
target_link_libraries( test22_PbarP7p3 ${Geant4_LIBRARIES} )

install( TARGETS test22_PbarP7p3 DESTINATION bin )

#-----
add_executable( test22_PbarP9p1 PbarP9p1.cc ${sources} ${headers} )
target_link_libraries( test22_PbarP9p1 ${Geant4_LIBRARIES} )

install( TARGETS test22_PbarP9p1 DESTINATION bin )
