# - CXX compile features for Intel to support C++14/17 on CMake from 3.8
#
#-----------------------------------------------------------------------
# Add compile features for Intel - should eventually be placed
# into a module, as it will need exporting for use by clients
if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  # CMake 3.8,3.9,3.10 support up to 16.0, c++14
  #       3.11,3,12,3.13,3.14,3.15,3.16 support up to 18, c++17
  # So need to provide addons to provide c++17 on 3.8-3.10
  if(CMAKE_VERSION VERSION_LESS 3.11)
    if("x${CMAKE_CXX_SIMULATE_ID}" STREQUAL "xMSVC")
      if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18.0.0)
        set(CMAKE_CXX17_STANDARD_COMPILE_OPTION "-Qstd=c++17")
        set(CMAKE_CXX17_EXTENSION_COMPILE_OPTION "-Qstd=c++17")

        list(APPEND CMAKE_CXX17_COMPILE_FEATURES cxx_std_17)
        set(CMAKE_CXX_COMPILE_FEATURES ${CMAKE_CXX_COMPILE_FEATURES} "${CMAKE_CXX17_COMPILE_FEATURES}")
      endif()
    else()
      if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18.0.0)
        set(CMAKE_CXX17_STANDARD_COMPILE_OPTION "-std=c++17")
        set(CMAKE_CXX17_EXTENSION_COMPILE_OPTION "-std=gnu++17")
        list(APPEND CMAKE_CXX17_COMPILE_FEATURES cxx_std_17)
        set(CMAKE_CXX_COMPILE_FEATURES ${CMAKE_CXX_COMPILE_FEATURES} "${CMAKE_CXX17_COMPILE_FEATURES}")
      endif()
    endif()
  endif()
endif()
