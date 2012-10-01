# Locate Pythia6 library
# in a directory defined via PYTHIA6 environment variable
#
# Defines:
#  PYTHIA6_FOUND
#  PYTHIA6_LIBRARIES

find_library(PYTHIA6_LIBRARY NAMES Pythia6 pythia6-$ENV{PYTHIA6_VERSION}
			     HINTS $ENV{PYTHIA6} $ENV{PYTHIA6}/lib)

set(PYTHIA6_LIBRARIES ${PYTHIA6_LIBRARY})                              
#message(STATUS PYTHIA6_LIBRARIES ${PYTHIA6_LIBRARIES} )

# handle the QUIETLY and REQUIRED arguments and set PYTHIA6_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Pythia6 DEFAULT_MSG PYTHIA6_LIBRARIES)

mark_as_advanced(PYTHIA6_FOUND PYTHIA6_LIBRARIES)
