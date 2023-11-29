if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|(Apple)?Clang")
    set(ptl_sanitize_args "-fsanitize=${PTL_SANITIZER_TYPE}")

    if(NOT ("${PTL_SANITIZER_TYPE}" STREQUAL "thread"))
        string(
            APPEND ptl_sanitize_args
            " -fno-optimize-sibling-calls -fno-omit-frame-pointer -fno-inline-functions")
    endif()
else()
    set(ptl_sanitize_args)
    message(
        WARNING
            "PTL_USE_SANITIZER is ON, but is not supported for compiler '${CMAKE_CXX_COMPILER_ID}"
        )
endif()

# We're setting the CXX flags and C flags beacuse they're propagated down independent of
# build type. Once CMake 3.13 is the minimum, we can use add_{compile,link}_options +
# genexs instead
string(APPEND CMAKE_CXX_FLAGS " ${ptl_sanitize_args}")
string(APPEND CMAKE_C_FLAGS " ${ptl_sanitize_args}")
string(APPEND CMAKE_EXE_LINKER_FLAGS " ${ptl_sanitize_args}")
string(APPEND CMAKE_SHARED_LINKER_FLAGS " ${ptl_sanitize_args}")
string(APPEND CMAKE_MODULE_LINKER_FLAGS " ${ptl_sanitize_args}")
