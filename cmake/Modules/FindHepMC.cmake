# - Locate HepMC library
# in a directory defined via  HEPMC_ROOT_DIR or HEPMC_DIR environment variable
# Defines:
#
#  HEPMC_FOUND
#  HEPMC_INCLUDE_DIR
#  HEPMC_INCLUDE_DIRS (not cached)
#  HEPMC_LIBRARIES
#  HEPMC_FIO_LIBRARIES

find_path(HEPMC_INCLUDE_DIR HepMC/GenEvent.h
          HINTS $ENV{HEPMC_ROOT_DIR}/include ${HEPMC_ROOT_DIR}/include
          $ENV{HEPMC_DIR}/include ${HEPMC_DIR}/include)

find_library(HEPMC_LIBRARIES NAMES HepMC
             HINTS $ENV{HEPMC_ROOT_DIR}/lib ${HEPMC_ROOT_DIR}/lib
             HINTS $ENV{HEPMC_DIR}/lib ${HEPMC_DIR}/lib)

get_filename_component(HEPMC_LIBRARY_DIR ${HEPMC_LIBRARIES} PATH)
set(HEPMC_FIO_LIBRARIES "-L${HEPMC_LIBRARY_DIR} -lHepMCfio")

set(HEPMC_INCLUDE_DIRS ${HEPMC_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments and set HEPMC_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(HepMC DEFAULT_MSG HEPMC_INCLUDE_DIR HEPMC_LIBRARIES)

mark_as_advanced(HEPMC_FOUND HEPMC_INCLUDE_DIR HEPMC_LIBRARIES)
