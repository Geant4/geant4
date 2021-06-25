#---Adding all geometry examples subdirectories explicitly

cmake_minimum_required(VERSION 3.12...3.20)

add_subdirectory(transforms)

find_package(VecGeom QUIET)
if(VecGeom_FOUND)
  add_subdirectory(vecGeomNavigation)
endif()
