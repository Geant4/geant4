#---Adding examples subdirectories explicitly
#   and a custom target to for building all basic examples -------------

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

add_subdirectory(B2b)

add_custom_target(tbb DEPENDS B2bTBB)
