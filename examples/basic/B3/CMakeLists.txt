#---Adding example B3 subdirectories explicitly
#   and a custom target to for building all example B2 options ----------
cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

add_subdirectory(B3a)
add_subdirectory(B3b)

add_custom_target(B3 DEPENDS exampleB3a exampleB3b)
