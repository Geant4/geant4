cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test46)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

if(GEANT4_ENABLE_TESTING)
  GEANT4_EXECUTABLE(test46 test46.cc src/*.cc)

  #---Test definitions-----------------------------------------------------
  GEANT4_ADD_TEST(test46 COMMAND test46 ${CMAKE_CURRENT_SOURCE_DIR}/test46.in
                       BUILD test46
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
  add_executable(test46 test46.cc ${sources} ${headers})
  target_link_libraries(test46 ${Geant4_LIBRARIES} )
  install(TARGETS test46 DESTINATION bin)

endif()
