cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test45)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

if(GEANT4_ENABLE_TESTING)
  GEANT4_EXECUTABLE(test45 test45.cc src/*.cc)

  #---Test definitions----------------------------------------------------------------------------------
  GEANT4_ADD_TEST(test45 COMMAND test45 ${CMAKE_CURRENT_SOURCE_DIR}/test45.in
                       BUILD test45
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
  add_executable(test45 test45.cc ${sources} ${headers})
  target_link_libraries(test45 ${Geant4_LIBRARIES} )
  install(TARGETS test45 DESTINATION bin)

endif()
