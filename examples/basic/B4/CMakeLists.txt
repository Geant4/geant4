
#---Adding example N3 subdirectories explicitly 
#   and a custom target to for building all example N3 options ----------

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

add_subdirectory(B4a)
add_subdirectory(B4b)
add_subdirectory(B4c)
add_subdirectory(B4d)

add_custom_target(B4 DEPENDS exampleB4a exampleB4b exampleB4c exampleB4d)
