#---Adding examples subdirectories explicitly
#   and a custom target to for building all basic examples -------------

cmake_minimum_required(VERSION 3.16...3.21)

add_subdirectory(B2b)

add_custom_target(tbb DEPENDS B2bTBB)
