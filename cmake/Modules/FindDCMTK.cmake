# Locate DCMTK package
# in a directory defined via DCMTK_DIR CMake or environment variable
# Defines:
#
#  DCMTK_INCLUDE_DIRS
#  DCMTK_LIBRARIES
#  DCMTK_FOUND

# includes
#
find_path(DCMTK_INCLUDE_DIR dcmtk/dcmdata/dcfilefo.h
          HINTS $ENV{DCMTK_DIR}/include ${DCMTK_DIR}/include)
set(DCMTK_INCLUDE_DIRS ${DCMTK_INCLUDE_DIR})


# libraries
#
find_library(DCMTK_DCMPSTAT_LIBRARY
             NAMES dcmpstat 
             HINTS "$ENV{DCMTK_DIR}/lib" "${DCMTK_DIR}/lib")

get_filename_component(DCMTK_LIBRARY_DIR 
	         ${DCMTK_DCMPSTAT_LIBRARY} DIRECTORY)

set(DCMTK_LIBRARIES 
	-L${DCMTK_LIBRARY_DIR} -ldcmpstat -ldcmwlm -lijg8 -ldcmdata -ldcmjpeg -ldcmqrdb 
	-li2d -loflog -ldcmdsig -ldcmjpls -ldcmsr -lijg12 -lofstd -ldcmimage -ldcmnet 
	-ldcmtls -lijg16 -ldcmjpeg -ldcmrt -lcharls -ldcmimgle)

# handle the QUIETLY and REQUIRED arguments and set DCMTK_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(DCMTK DEFAULT_MSG DCMTK_INCLUDE_DIRS DCMTK_LIBRARIES)

mark_as_advanced(DCMTK_FOUND DCMTK_INCLUDE_DIRS DCMTK_LIBRARIES)
