cmake_minimum_required( VERSION 2.6 FATAL_ERROR )
project( test22_hA )
find_package( Geant4 REQUIRED )

include( ${Geant4_USE_FILE} )
include_directories( ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/../src ${PROJECT_SOURCE_DIR}/../include)
file( GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/../src/*.cc)
file( GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/../include/*.hh)

#-----
add_executable( test22_E941 E941.cc ${sources} ${headers} )
target_link_libraries( test22_E941 ${Geant4_LIBRARIES} )

install( TARGETS test22_E941 DESTINATION bin )

#-----
add_executable( test22_hA100 hA100.cc ${sources} ${headers} )
target_link_libraries( test22_hA100 ${Geant4_LIBRARIES} )

install( TARGETS test22_hA100 DESTINATION bin )

#-----
add_executable( test22_hA120 hA120.cc ${sources} ${headers} )
target_link_libraries( test22_hA120 ${Geant4_LIBRARIES} )

install( TARGETS test22_hA120 DESTINATION bin )

#-----
add_executable( test22_HARP HARP.cc ${sources} ${headers} )
target_link_libraries( test22_HARP ${Geant4_LIBRARIES} )

install( TARGETS test22_HARP DESTINATION bin )

#-----
add_executable( test22_NA35 NA35.cc ${sources} ${headers} )
target_link_libraries( test22_NA35 ${Geant4_LIBRARIES} )

install( TARGETS test22_NA35 DESTINATION bin )

#-----
add_executable( test22_NA49 NA49.cc ${sources} ${headers} )
target_link_libraries( test22_NA49 ${Geant4_LIBRARIES} )

install( TARGETS test22_NA49 DESTINATION bin )

#-----
add_executable( test22_NA61 NA61.cc ${sources} ${headers} )
target_link_libraries( test22_NA61 ${Geant4_LIBRARIES} )

install( TARGETS test22_NA61 DESTINATION bin )

#-----
add_executable( test22_pA200 pA200.cc ${sources} ${headers} )
target_link_libraries( test22_pA200 ${Geant4_LIBRARIES} )

install( TARGETS test22_pA200 DESTINATION bin )

#---Test definition-----------------------------------------

GEANT4_ADD_TEST( test22_NA49-build BUILD test22_NA49 
                       LABELS Nightly Continuous ) 

GEANT4_ADD_TEST( test22_NA49-exerciser
                       COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test22_NA49 ${CMAKE_CURRENT_SOURCE_DIR}/NA49pC.mac
                       DEPENDS test22_NA49-build
                       LABELS Nightly Continuous
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} )
