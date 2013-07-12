cmake_minimum_required( VERSION 2.6 FATAL_ERROR )
project( test22_PPinclusive )
find_package( Geant4 REQUIRED )

include( ${Geant4_USE_FILE} )
include_directories( ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/../src ${PROJECT_SOURCE_DIR}/../include)
file( GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/../src/*.cc)
file( GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/../include/*.hh)

#-----
add_executable( test22_PP100 PP100.cc ${sources} ${headers} )
target_link_libraries( test22_PP100 ${Geant4_LIBRARIES} )

install( TARGETS test22_PP100 DESTINATION bin )

#-----
add_executable( test22_PP12 PP12.cc ${sources} ${headers} )
target_link_libraries( test22_PP12 ${Geant4_LIBRARIES} )

install( TARGETS test22_PP12 DESTINATION bin )

#-----
add_executable( test22_PP158 PP158.cc ${sources} ${headers} )
target_link_libraries( test22_PP158 ${Geant4_LIBRARIES} )

install( TARGETS test22_PP158 DESTINATION bin )

#-----
add_executable( test22_PP175 PP175.cc ${sources} ${headers} )
target_link_libraries( test22_PP175 ${Geant4_LIBRARIES} )

install( TARGETS test22_PP175 DESTINATION bin )

#-----
add_executable( test22_PP205 PP205.cc ${sources} ${headers} )
target_link_libraries( test22_PP205 ${Geant4_LIBRARIES} )

install( TARGETS test22_PP205 DESTINATION bin )

#-----
add_executable( test22_PP24 PP24.cc ${sources} ${headers} )
target_link_libraries( test22_PP24 ${Geant4_LIBRARIES} )

install( TARGETS test22_PP24 DESTINATION bin )

#-----
add_executable( test22_PP360 PP360.cc ${sources} ${headers} )
target_link_libraries( test22_PP360 ${Geant4_LIBRARIES} )

install( TARGETS test22_PP360 DESTINATION bin )

#-----
add_executable( test22_PP400 PP400.cc ${sources} ${headers} )
target_link_libraries( test22_PP400 ${Geant4_LIBRARIES} )

install( TARGETS test22_PP400 DESTINATION bin )

#-----
add_executable( test22_PP69 PP69.cc ${sources} ${headers} )
target_link_libraries( test22_PP69 ${Geant4_LIBRARIES} )

install( TARGETS test22_PP69 DESTINATION bin )
