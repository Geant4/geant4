# Locate CERNLIB libraries (needed for g4tools/hbook)
# in a directory defined via CERNLIB environment variable
#
# Defines:
#  HBOOK_FOUND
#  HBOOK_LIBRARIES

find_library(PACKLIB_LIBRARY NAMES packlib
             HINTS $ENV{CERNLIB} $ENV{CERNLIB}/lib)

# handle the QUIETLY and REQUIRED arguments and set HBOOK_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(HBOOK DEFAULT_MSG PACKLIB_LIBRARY)

if (HBOOK_FOUND)
  get_filename_component(CERNLIB_LIBRARY_DIR ${PACKLIB_LIBRARY} PATH)
  set(HBOOK_LIBRARIES "-L${CERNLIB_LIBRARY_DIR} -lpacklib -lmathlib -lgfortran -lcrypt")
endif()    

mark_as_advanced(HBOOK_FOUND HBOOK_LIBRARIES)
