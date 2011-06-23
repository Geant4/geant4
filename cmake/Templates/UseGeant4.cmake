# - Use file for Geant4
# This file should be included after a find_package call has successfully
# located Geant4.
#
#  Geant4_USE_FILE : Point to the location of the use file for the found
#                    Geant4 installation.
#
# Inclusion of this file results in the addition of the Geant4 compile
# definitions and include directories to those of the directory in which
# this files is included.
#
# The recommended Geant4 compiler flags are also added to CMAKE_CXX_FLAGS
# Advanced users requiring special sets of flags should therefore not use this
# file, preferring the direct use of the Geant4_XXXX variables set by the
# Geant4Config file.
#
# The use file also defines a simple macro to help with collating sources
# for simple user applications.
#
#  macro GEANT4_COLLATE_APPLICATION_SOURCES(sources)
#        Create the list of sources in a standard Geant4 application example.
#        If the project is layed out as
#        +- project
#           +- CMakeLists.txt
#           +- project.cc
#           +- include/
#           |  +- project_impl_a.hh
#           |  +- ...
#           +- src/
#              +- project_impl_a.cc
#              +- ...
#        Then when called in the CMakeLists.txt under 'project', it will
#        collate all the .hh and .cc files into the sources variable, and
#        add the include directory to those to be searched for.
#        The CMakeLists.txt for project would then read
#       
#        cmake_minimum_required(VERSION 2.6.4)
#        project(myproject)
#        find_package(Geant4 REQUIRED)
#        include(${Geant4_USE_FILE})
#        GEANT4_COLLATE_SOURCES(my_project_sources)
#        add_executable(myproject ${my_project_sources})
#        target_link_libraries(myproject ${Geant4_LIBRARIES})
#
#        This gives you maximum flexibility in setting up your project.
#

#----------------------------------------------------------------------------
# macro GEANT4_COLLATE_SOURCES(source_dest_var)
#
macro(GEANT4_COLLATE_APPLICATION_SOURCES source_dest_var)
    file(GLOB_RECURSE 
        ${source_dest_var} 
        ${CMAKE_CURRENT_SOURCE_DIR}/*.hh 
        ${CMAKE_CURRENT_SOURCE_DIR}*.cc
    )
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
endmacro()


#----------------------------------------------------------------------------
# We need to set the compile definitions and include directories
#
add_definitions(${Geant4_DEFINITIONS})
include_directories(${Geant4_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Because Geant4 is sensitive to the compiler flags, lets set the base set
# here. This reproduces as far as possible the behaviour of the original
# makefile system.
#
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Geant4_CXX_FLAGS}")


