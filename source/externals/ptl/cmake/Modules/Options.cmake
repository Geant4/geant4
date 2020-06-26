
################################################################################
#
#        PTL Options
#
################################################################################

include(MacroUtilities)

set(_FEATURE )
if(NOT PTL_MASTER_PROJECT)
    set(_FEATURE NO_FEATURE)
endif()

set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Require the C++ standard" FORCE)
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Disable GNU extensions")

# features
add_feature(CMAKE_BUILD_TYPE "Build type (Debug, Release, RelWithDebInfo, MinSizeRel)")
add_feature(CMAKE_INSTALL_PREFIX "Installation prefix")
add_feature(CMAKE_CXX_STANDARD "C++11 STL standard")
add_feature(${PROJECT_NAME}_C_FLAGS "C flags for project")
add_feature(${PROJECT_NAME}_CXX_FLAGS "C++ flags for project")

# options (always available)
add_option(BUILD_STATIC_LIBS "Build static library" ON ${_FEATURE})
add_option(BUILD_SHARED_LIBS "Build shared library" ON ${_FEATURE})
add_option(PTL_BUILD_EXAMPLES "Build examples" OFF ${_FEATURE})
add_option(PTL_BUILD_DOCS "Build documentation with Doxygen" OFF ${_FEATURE})
add_option(PTL_DEVELOPER_INSTALL "Install headers, cmake export, and shared libs" ON ${_FEATURE})

add_option(PTL_USE_TBB "Enable TBB" ON ${_FEATURE})
add_option(PTL_USE_GPU "Enable GPU preprocessor" OFF ${_FEATURE})
add_option(PTL_USE_SANITIZER "Enable -fsanitize=<type>" OFF ${_FEATURE})
add_option(PTL_USE_CLANG_TIDY "Enable running clang-tidy on" OFF ${_FEATURE})
add_option(PTL_USE_COVERAGE "Enable code coverage" OFF ${_FEATURE})
add_option(PTL_USE_PROFILE "Enable profiling" OFF ${_FEATURE})

if(PTL_USE_ARCH)
    add_option(PTL_USE_AVX512 "Enable AVX-512 flags (if available)" OFF ${_FEATURE})
endif()

if(PTL_USE_SANITIZER)
    add_feature(PTL_SANITIZER_TYPE "Sanitizer type (-fsanitize=<type>)")
    set(PTL_SANITIZER_TYPE leak CACHE STRING "Sanitizer type (-fsanitize=<type>)")
endif()

if(PTL_USE_GPU)
    # Check if CUDA can be enabled
    find_package(CUDA QUIET)
    if(CUDA_FOUND)
        check_language(CUDA)
        if(CMAKE_CUDA_COMPILER)
            enable_language(CUDA)
        else()
            message(STATUS "No CUDA support")
            set(_USE_CUDA OFF)
        endif()
    else()
        set(_USE_CUDA OFF)
        set(PTL_USE_GPU OFF)
    endif()

    if(_USE_CUDA)
        list(APPEND ${PROJECT_NAME}_DEFINITIONS PTL_USE_GPU)
        add_option(PTL_USE_CUDA "Enable CUDA option for GPU execution" ${_USE_CUDA} ${_FEATURE})
        add_option(PTL_USE_NVTX "Enable NVTX for Nsight" ${_USE_CUDA} ${_FEATURE})

        set(CUDA_ARCH "sm_35" CACHE STRING "CUDA architecture flag")
        if(_FEATURE)
            add_feature(CUDA_ARCH "CUDA architecture (e.g. sm_35)")
            add_feature(CUDA_GENERATED_OUTPUT_DIR "CUDA output directory for generated files")
        endif()
    endif()

endif()

# RPATH settings
set(_RPATH_LINK OFF)
if(APPLE)
    set(_RPATH_LINK ON)
endif()
add_option(CMAKE_INSTALL_RPATH_USE_LINK_PATH "Hardcode installation rpath based on link path" ${_RPATH_LINK} ${_FEATURE})
unset(_RPATH_LINK)

# clang-tidy
if(PTL_USE_CLANG_TIDY)
    find_program(CLANG_TIDY_COMMAND NAMES clang-tidy)
    add_feature(CLANG_TIDY_COMMAND "Path to clang-tidy command")
    if(NOT CLANG_TIDY_COMMAND)
        message(WARNING "PTL_USE_CLANG_TIDY is ON but clang-tidy is not found!")
        set(PTL_USE_CLANG_TIDY OFF)
    else()
        set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND}")

        # Create a preprocessor definition that depends on .clang-tidy content so
        # the compile command will change when .clang-tidy changes.  This ensures
        # that a subsequent build re-runs clang-tidy on all sources even if they
        # do not otherwise need to be recompiled.  Nothing actually uses this
        # definition.  We add it to targets on which we run clang-tidy just to
        # get the build dependency on the .clang-tidy file.
        file(SHA1 ${PROJECT_SOURCE_DIR}/.clang-tidy clang_tidy_sha1)
        set(CLANG_TIDY_DEFINITIONS "CLANG_TIDY_SHA1=${clang_tidy_sha1}")
        unset(clang_tidy_sha1)
    endif()
    configure_file(${PROJECT_SOURCE_DIR}/.clang-tidy ${PROJECT_SOURCE_DIR}/.clang-tidy COPYONLY)
endif()
