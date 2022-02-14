
#------------------------------------------------------------------------------#
#   CUDA function
#
function(NVCUDA_COMPILE_PTX)
    set(options "")
    set(oneValueArgs TARGET_PATH GENERATED_FILES OUTPUT_DIR)
    set(multiValueArgs NVCC_OPTIONS SOURCES)
    cmake_parse_arguments(NVCUDA_COMPILE_PTX "${options}" "${oneValueArgs}"
        "${multiValueArgs}" ${ARGN})

    # Match the bitness of the ptx to the bitness of the application
    set( MACHINE "--machine=32" )
    if( CMAKE_SIZEOF_VOID_P EQUAL 8)
        set( MACHINE "--machine=64" )
    endif()

    set(_OUTPUT_DIR ${NVCUDA_COMPILE_PTX_OUTPUT_DIR})
    if("${_OUTPUT_DIR}" STREQUAL "")
        set(_OUTPUT_DIR ${CUDA_GENERATED_OUTPUT_DIR})
    endif()

    # Custom build rule to generate ptx files from cuda files
    foreach(input ${NVCUDA_COMPILE_PTX_SOURCES})
        get_filename_component(input_we ${input} NAME_WE)

        # generate the *.ptx files inside "ptx" folder inside
        # the executable's output directory.
        set( output "${_OUTPUT_DIR}/${input_we}.ptx" )

        list(APPEND PTX_FILES  ${output})

        add_custom_command(
            OUTPUT  ${output}
            DEPENDS ${input}
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            COMMAND ${CUDA_NVCC_EXECUTABLE} ${MACHINE}
                --ptx ${NVCUDA_COMPILE_PTX_NVCC_OPTIONS} ${input}
                -o ${output}
        )
    endforeach()

    set(${NVCUDA_COMPILE_PTX_GENERATED_FILES} ${PTX_FILES} PARENT_SCOPE)
endfunction()
