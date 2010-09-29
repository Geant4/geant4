# - Geant4MakeRules_fortran
# Sets the default make rules for a Fortran build, specifically the 
# initialization of the compiler flags on a platform dependent basis
#

message(STATUS "setting default compiler flags for Fortran")

# We have to rely on the compiler name
get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

# gfortran...
if(Fortran_COMPILER_NAME STREQUAL "gfortran")
    set(CMAKE_Fortran_FLAGS_INIT "-fno-automatic -fno-backslash -fno-second-underscore")
    set(CMAKE_Fortran_FLAGS_DEBUG_INIT "-g") 
    set(CMAKE_Fortran_FLAGS_RELEASE_INIT "-O3")
    set(CMAKE_Fortran_FLAGS_MINSIZEREL_INIT "-Os")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")
elseif(Fortran_COMPILER_NAME STREQUAL "ifort")
    set(CMAKE_Fortran_FLAGS_INIT "-noautomatic -assume nobscc -assume no2underscores")
    set(CMAKE_Fortran_FLAGS_DEBUG_INIT "-g") 
    set(CMAKE_Fortran_FLAGS_RELEASE_INIT "-O3")
    set(CMAKE_Fortran_FLAGS_MINSIZEREL_INIT "-Os")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")
endif()

# ifort...
