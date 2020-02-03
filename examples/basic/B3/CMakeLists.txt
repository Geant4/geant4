
#---Adding example B2 subdirectories explicitly 
#   and a custom target to for building all example B2 options ----------

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

add_subdirectory(B3a)
add_subdirectory(B3b)

add_custom_target(B3 DEPENDS exampleB3a exampleB3b)
