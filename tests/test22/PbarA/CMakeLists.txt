cmake_minimum_required( VERSION 2.6 FATAL_ERROR )
project( test22_PbarA )
find_package( Geant4 REQUIRED )

include( ${Geant4_USE_FILE} )
include_directories( ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/../src ${PROJECT_SOURCE_DIR}/../include)
file( GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/../src/*.cc)
file( GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/../include/*.hh)

#-----
add_executable( test22_Nb_C_750 Nb_C_750.cc ${sources} ${headers} )
target_link_libraries( test22_Nb_C_750 ${Geant4_LIBRARIES} )

install( TARGETS test22_Nb_C_750 DESTINATION bin )

#-----
add_executable( test22_PbA_1p22_nE PbA_1p22_nE.cc ${sources} ${headers} )
target_link_libraries( test22_PbA_1p22_nE ${Geant4_LIBRARIES} )

install( TARGETS test22_PbA_1p22_nE DESTINATION bin )

#-----
add_executable( test22_PbA_Rest PbA_Rest.cc ${sources} ${headers} )
target_link_libraries( test22_PbA_Rest ${Geant4_LIBRARIES} )

install( TARGETS test22_PbA_Rest DESTINATION bin )

#-----
add_executable( test22_PbC_608 PbC_608.cc ${sources} ${headers} )
target_link_libraries( test22_PbC_608 ${Geant4_LIBRARIES} )

install( TARGETS test22_PbC_608 DESTINATION bin )

#-----
add_executable( test22_PbDNe_Rest PbDNe_Rest.cc ${sources} ${headers} )
target_link_libraries( test22_PbDNe_Rest ${Geant4_LIBRARIES} )

install( TARGETS test22_PbDNe_Rest DESTINATION bin )

#-----
add_executable( test22_PbTa_4 PbTa_4.cc ${sources} ${headers} )
target_link_libraries( test22_PbTa_4 ${Geant4_LIBRARIES} )

install( TARGETS test22_PbTa_4 DESTINATION bin )

#-----
add_executable( test22_PbU_608 PbU_608.cc ${sources} ${headers} )
target_link_libraries( test22_PbU_608 ${Geant4_LIBRARIES} )

install( TARGETS test22_PbU_608 DESTINATION bin )
