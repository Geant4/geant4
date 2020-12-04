#---Adding all runAndEvent examples subdirectories explicitly 

cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

add_subdirectory(RE01)
add_subdirectory(RE02)
add_subdirectory(RE03)
add_subdirectory(RE04)
add_subdirectory(RE05)
add_subdirectory(RE06)
