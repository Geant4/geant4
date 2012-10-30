#---Adding examples subdirectories explicitly 
#   and a custom target to for building all basic examples -------------

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

add_subdirectory(N01)
add_subdirectory(N02)
add_subdirectory(N03)
add_subdirectory(N04)
add_subdirectory(N05)
add_subdirectory(N06)
add_subdirectory(N07)

add_custom_target(novice DEPENDS N01 N02 N03 N04 N05 N06 N07)
