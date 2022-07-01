# Only search for clang-tidy if user has not supplied it through CMAKE_CXX_CLANG_TIDY
if(NOT DEFINED CMAKE_CXX_CLANG_TIDY)
    find_program(
        PTL_CLANG_TIDY_COMMAND NAMES clang-tidy clang-tidy-12 clang-tidy-11 clang-tidy-10
                                     clang-tidy-9 clang-tidy-8 clang-tidy-7)
    mark_as_advanced(PTL_CLANG_TIDY_COMMAND)
    if(NOT PTL_CLANG_TIDY_COMMAND)
        message(
            FATAL_ERROR
                "PTL_USE_CLANG_TIDY is ON, but no clang-tidy program could be found")
    endif()
    set(CMAKE_CXX_CLANG_TIDY "${PTL_CLANG_TIDY_COMMAND}")
endif()

# Create a preprocessor definition that depends on .clang-tidy content so the compile
# command will change when .clang-tidy changes. Nothing actually uses this definition, it
# is simply added to trigger a rebuild when *only* .clang-tidy changes.
file(SHA1 ${PROJECT_SOURCE_DIR}/.clang-tidy __ptl_clang_tidy_sha1)
add_compile_definitions("PTL_CLANG_TIDY_SHA1=${__ptl_clang_tidy_sha1}")
unset(__ptl_clang_tidy_sha1)
configure_file(${PROJECT_SOURCE_DIR}/.clang-tidy ${PROJECT_SOURCE_DIR}/.clang-tidy
               COPYONLY)
