
#---Adding example B2 subdirectories explicitly 
#   and a custom target to for building all example B2 options ----------

cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

add_subdirectory(gflash1)
add_subdirectory(gflash2)
add_subdirectory(gflash3)
add_subdirectory(gflasha)

add_custom_target(gflash DEPENDS gflash1 gflash2 gflash3 gflasha)
