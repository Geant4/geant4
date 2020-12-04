#---Adding all optical examples subdirectories explicitly 

cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

add_subdirectory(OpNovice)
add_subdirectory(OpNovice2)
add_subdirectory(LXe)
add_subdirectory(wls)
