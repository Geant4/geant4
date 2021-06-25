# Locate Pythia8 library
# in a directory defined via PYTHIA8 environment variable
#
# Defines:
#  PYTHIA8_FOUND
#  PYTHIA8_LIBRARIES
#  PYTHIA8_INCLUDES

find_library(PYTHIA8_LIBRARY NAMES pythia8
			     HINTS $ENV{PYTHIA8} $ENV{PYTHIA8}/lib)

set(PYTHIA8_LIBRARIES ${PYTHIA8_LIBRARY})
#message(STATUS PYTHIA8_LIBRARIES ${PYTHIA8_LIBRARIES} )

find_path( PYTHIA8_INCLUDES Pythia8/Pythia.h
           HINTS $ENV{PYTHIA8} $ENV{PYTHIA8}/include $ENV{PYTHIA8}/include/Pythia8 )
set(PYTHIA8_INCLUDES ${PYTHIA8_INCLUDES})

# handle the QUIETLY and REQUIRED arguments and set PYTHIA8_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Pythia8 DEFAULT_MSG PYTHIA8_LIBRARIES)

mark_as_advanced(PYTHIA8_FOUND PYTHIA8_LIBRARIES PYTHIA8_INCLUDES)
