#---Adding all biasing examples subdirectories explicitly 

cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

add_subdirectory(B01)
add_subdirectory(B02)
add_subdirectory(B03)
add_subdirectory(GB01)
add_subdirectory(GB02)
add_subdirectory(GB03)
add_subdirectory(GB04)
add_subdirectory(GB05)
add_subdirectory(GB06)
add_subdirectory(GB07)
add_subdirectory(ReverseMC01)
