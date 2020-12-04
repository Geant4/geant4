#---Adding all hadronic examples subdirectories explicitly 

cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

add_subdirectory(Hadr00)
add_subdirectory(Hadr01)
add_subdirectory(Hadr02)
add_subdirectory(Hadr03)
add_subdirectory(Hadr04)
add_subdirectory(Hadr06)
add_subdirectory(Hadr07)
add_subdirectory(Hadr08)
add_subdirectory(Hadr09)
add_subdirectory(FissionFragment)
add_subdirectory(NeutronSource)
