# - Locate LHAPDF library
# in a directory defined via  LHAPDF_ROOT_DIR or LHAPDF_DIR environment variable
# Defines:
#
#  LHAPDF_FOUND
#  LHAPDF_INCLUDE_DIR
#  LHAPDF_INCLUDE_DIRS (not cached)
#  LHAPDF_LIBRARIES

find_path(LHAPDF_INCLUDE_DIR LHAPDF.h
          HINTS $ENV{LHAPDF_ROOT_DIR}/include ${LHAPDF_ROOT_DIR}/include
          $ENV{LHAPDF_DIR}/include ${LHAPDF_DIR}/include)

find_library(LHAPDF_LIBRARIES NAMES LHAPDF
             HINTS $ENV{LHAPDF_ROOT_DIR}/lib ${LHAPDF_ROOT_DIR}/lib
             HINTS $ENV{LHAPDF_DIR}/lib ${LHAPDF_DIR}/lib
             HINTS $ENV{LHAPDF_ROOT_DIR}/lib64 ${LHAPDF_ROOT_DIR}/lib64
             HINTS $ENV{LHAPDF_DIR}/lib64 ${LHAPDF_DIR}/lib64
             )


set(LHAPDF_INCLUDE_DIRS ${LHAPDF_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments and set LHAPDF_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LHAPDF DEFAULT_MSG LHAPDF_INCLUDE_DIR LHAPDF_LIBRARIES)

mark_as_advanced(LHAPDF_FOUND LHAPDF_INCLUDE_DIR LHAPDF_LIBRARIES)
