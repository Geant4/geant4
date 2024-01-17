#---Adding all geometry examples subdirectories explicitly

cmake_minimum_required(VERSION 3.16...3.27)

add_subdirectory(transforms)

find_package(VecGeom QUIET)
if(VecGeom_FOUND)
  add_subdirectory(vecGeomNavigation)
endif()
