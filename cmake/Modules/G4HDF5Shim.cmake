# - G4HDF5Shim
#
# Geant4's Geant4Config.cmake file aims to support CMake 3.8 and newer
# The HDF5 dependency is located through CMake's builtin FindHDF5
# module, but this does not supply imported targets until CMake 3.20.
# It may use HDF5's hdf5-config.cmake file if available, so create
# custom imported target hdf5::hdf5 to allow both cases to be handled
# without interference with either.
 
if(HDF5_FOUND)
 # If we're in MT mode, found HDF5 must also support MT
 if(GEANT4_BUILD_MULTITHREADED OR Geant4_multithreaded_FOUND)
   include(CheckCXXSymbolExists)
   set(CMAKE_REQUIRED_INCLUDES "${HDF5_INCLUDE_DIRS}")
   check_cxx_symbol_exists(H5_HAVE_THREADSAFE "H5pubconf.h" GEANT4_HAVE_H5_HAVE_THREADSAFE)
   unset(CMAKE_REQUIRED_INCLUDES)

   if(NOT GEANT4_HAVE_H5_HAVE_THREADSAFE)
     message(FATAL_ERROR
       "Found an install of HDF5, but it was not built with support for thread safety. "
       "Either build Geant4 in single threaded mode, or use/reinstall HDF5 with "
       "thread safety enabled. See HDF5's install guides, available from https://support.hdfgroup.org/HDF5/release/, for instructions on this.\n"
       )
   endif()
 endif()

 # If FindHDF5 does not yet supply imported targets, we
 # create an internal INTERFACE target to wrap these.
 # This still hard-codes include/library paths, but limits it
 # to one place. Later, we'll create proper imported targets
 # with re-finds but for now this is the best minimally invasive proceedure
 if(NOT TARGET hdf5::hdf5)
   add_library(hdf5::hdf5 IMPORTED UNKNOWN)
   set_target_properties(hdf5::hdf5 PROPERTIES
     IMPORTED_LINK_INTERFACE_LANGUAGES "C"
     IMPORTED_LOCATION "${HDF5_C_LIBRARY_hdf5}"
     INTERFACE_INCLUDE_DIRECTORIES "${HDF5_C_INCLUDE_DIRS}"
     INTERFACE_LINK_LIBRARIES "${HDF5_C_LIBRARIES}"
     )
 endif()
endif()
