cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test30)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

if(GEANT4_ENABLE_TESTING)
  GEANT4_EXECUTABLE(test30 test30.cc src/*.cc)

  #---Test definitions----------------------------------------------------------------------------------
  GEANT4_ADD_TEST(test30 COMMAND test30 ${CMAKE_CURRENT_SOURCE_DIR}/test30.in
                       BUILD test30
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
else()
  # allow a build without testing enabled....
  # Locate sources and headers for this project
  include_directories(${PROJECT_SOURCE_DIR}/include 
                  	 ${Geant4_INCLUDE_DIR})
  file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
  file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

  #----------------------------------------------------------------------------
  # Add the executable, and link it to the Geant4 libraries
  #
  add_executable(test30 test30.cc ${sources} ${headers})
  target_link_libraries(test30 ${Geant4_LIBRARIES} )
  install(TARGETS test30 DESTINATION bin)

endif()
