#.rst:
# MSVC Compile Features for Geant4
# --------------------------------
#
# Due to use of C++17 features in visualization on Windows, Geant4 requires
# use of VS 2017 update 3 or newer. CMake 3.10 and higher support this out
# of the box, but with our minimum requirement of 3.8, we add this shim
# to:
#
# 1) Check the version requirement
# 2) Set the compile flags/features when using CMake < 3.10
#
if(MSVC)
  # This seems to indicate CMAKE_C_SIMULATE_ID
  if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  # Require MSVC that supports standard flags and std::filesystem
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 19.20)
      message(FATAL_ERROR "Geant4 requires MSVC 19.20 (Visual Studio 2019 Version 16.0) or newer")
    else()
      # Set cxxstd flags on CMake < 3.10
      if(CMAKE_VERSION VERSION_LESS 3.10)
        # VS 2015 Update 3 and above support language standard level flags,
        # with the default and minimum level being C++14.
        set(CMAKE_CXX98_STANDARD_COMPILE_OPTION "")
        set(CMAKE_CXX98_EXTENSION_COMPILE_OPTION "")
        set(CMAKE_CXX11_STANDARD_COMPILE_OPTION "")
        set(CMAKE_CXX11_EXTENSION_COMPILE_OPTION "")
        set(CMAKE_CXX14_STANDARD_COMPILE_OPTION "-std:c++14")
        set(CMAKE_CXX14_EXTENSION_COMPILE_OPTION "-std:c++14")
        set(CMAKE_CXX17_STANDARD_COMPILE_OPTION "-std:c++17")
        set(CMAKE_CXX17_EXTENSION_COMPILE_OPTION "-std:c++17")
        list(APPEND CMAKE_CXX17_COMPILE_FEATURES cxx_std_17)
        set(CMAKE_CXX_COMPILE_FEATURES ${CMAKE_CXX_COMPILE_FEATURES} "${CMAKE_CXX17_COMPILE_FEATURES}")
      endif()
    endif()
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  else()
    message(FATAL_ERROR "Geant4 requires Visual Studio or Clang.")
  endif()
endif()

