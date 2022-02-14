#---Adding example B3 subdirectories explicitly
#   and a custom target to for building all example B2 options ----------
cmake_minimum_required(VERSION 3.16...3.21)

add_subdirectory(B3a)
add_subdirectory(B3b)

add_custom_target(B3 DEPENDS exampleB3a exampleB3b)
