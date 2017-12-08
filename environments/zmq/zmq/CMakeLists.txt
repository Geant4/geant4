# - build G4zeromq library

include_directories(${PROJECT_SOURCE_DIR} ${GEANT4_INCLUDE_DIR}
                    ${ZMQ_INCLUDE_DIRS})

# library
add_library(
  G4zeromq SHARED
  G4ZMQServer.cc
)
# headers
set(HEADERS
  G4ZMQServer.hh
)

# install
install(TARGETS G4zeromq LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
