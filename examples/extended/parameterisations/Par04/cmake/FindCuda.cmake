# Find the CUDA include directory and library.
#
# This module defines the `cuda` imported target that encodes all
# necessary information in its target properties.
#
# This package is necessary for GPU memory profiling

find_library(
    Cuda_LIB
    NAMES cudart
    PATH_SUFFIXES lib lib32 lib64
    DOC "Cuda Runtime Library required for GPU Memory usage info"
)

find_path(
    Cuda_INCLUDE
    NAMES cuda.h
    PATH_SUFFIXES include
    DOC "Cuda Include directory required for GPU Memory usage info"
)

find_path(
    Cuda_Runtime_INCLUDE
    NAMES cuda_runtime_api.h
    PATH_SUFFIXES include
    DOC "Cuda Runtime Include directory required for GPU Memory usage info"
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    Cuda
    REQUIRED_VARS Cuda_LIB Cuda_INCLUDE Cuda_Runtime_INCLUDE
)

add_library(Cuda SHARED IMPORTED)
set_property(TARGET Cuda PROPERTY IMPORTED_LOCATION ${Cuda_LIB})
set_property(TARGET Cuda PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${Cuda_INCLUDE})
set_property(TARGET Cuda PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${Cuda_Runtime_INCLUDE})

mark_as_advanced(Cuda_FOUND Cuda_LIB Cuda_INCLUDE Cuda_Runtime_INCLUDE)
