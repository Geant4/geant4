if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|(Apple)?Clang")
    set(ptl_coverage_args "--coverage")
    check_cxx_compiler_flag("-fprofile-abs-path ${PTL_COVERAGE_FLAGS}"
                            PTL_HAS_PROFILE_ABS_PATH)
    if(PTL_HAS_PROFILE_ABS_PATH)
        string(APPEND ptl_coverage_args " -fprofile-abs-path")
    endif()
else()
    set(ptl_coverage_args)
    message(
        WARNING
            "PTL_USE_COVERAGE is ON, but is not supported for compiler '${CMAKE_CXX_COMPILER_ID}"
        )
endif()

# We're setting the DEBUG CXX flags and C flags because we only use coverage in Debug
# builds but want them to propagate down Once CMake 3.13 is the minimum, we can use
# add_{compile,link}_options + genexs instead
string(APPEND CMAKE_C_FLAGS_DEBUG " ${ptl_coverage_args}")
string(APPEND CMAKE_CXX_FLAGS_DEBUG " ${ptl_coverage_args}")
string(APPEND CMAKE_EXE_LINKER_FLAGS_DEBUG " ${ptl_coverage_args}")
string(APPEND CMAKE_SHARED_LINKER_FLAGS_DEBUG " ${ptl_coverage_args}")
string(APPEND CMAKE_MODULE_LINKER_FLAGS_DEBUG " ${ptl_coverage_args}")
