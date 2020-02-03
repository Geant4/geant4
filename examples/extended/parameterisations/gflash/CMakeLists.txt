
#---Adding example B2 subdirectories explicitly 
#   and a custom target to for building all example B2 options ----------

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

add_subdirectory(gflash1)
add_subdirectory(gflash2)
add_subdirectory(gflash3)
add_subdirectory(gflasha)

add_custom_target(gflash DEPENDS gflash1 gflash2 gflash3 gflasha)
