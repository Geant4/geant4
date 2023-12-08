#---Adding example B2 subdirectories explicitly
#   and a custom target to for building all example B2 options ----------
cmake_minimum_required(VERSION 3.16...3.27)

add_subdirectory(B2a)
add_subdirectory(B2b)

add_custom_target(B2 DEPENDS exampleB2a exampleB2b)
