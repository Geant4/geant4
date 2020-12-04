#---Adding example N3 subdirectories explicitly
#   and a custom target to for building all example N3 options ----------
cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

add_subdirectory(B4a)
add_subdirectory(B4b)
add_subdirectory(B4c)
add_subdirectory(B4d)

add_custom_target(B4 DEPENDS exampleB4a exampleB4b exampleB4c exampleB4d)
