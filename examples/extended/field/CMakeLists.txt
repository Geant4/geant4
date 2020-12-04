#---Adding all field examples subdirectories explicitly 

cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

add_subdirectory(BlineTracer)
add_subdirectory(field01)
add_subdirectory(field02)
add_subdirectory(field03)
add_subdirectory(field04)
add_subdirectory(field05)
add_subdirectory(field06)
