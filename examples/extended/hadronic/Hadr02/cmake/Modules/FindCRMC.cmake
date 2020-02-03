# - Locate CRMC library
# in a directory defined via  CRMCROOT environment variable
# Defines:  CRMC_FOUND
#           CRMC_INCLUDE_DIR
#           CRMC_LIBRARIES 
#           CRMC_LIBRARY_DIR
#           CRMC_BUILD_DIR
# Note: this file could eventually be moved into geant4/cmake/Modules/ .

# look if an environment variable CRMCROOT exists
# -DCRMCROOT
set(CRMCROOT $ENV{CRMCROOT})

find_library(CRMC_LIBRARIES libCrmc.so PATHS ${CRMCROOT}/Build/lib)
if (CRMC_LIBRARIES)
   set(CRMC_FOUND TRUE)
   set(CRMC_INCLUDE_DIR ${CRMCROOT}/src ${CRMCROOT}/Build/src)
   set(CRMC_LIBRARY_DIR ${CRMCROOT}/Build/lib)
   set(CRMC_BUILD_DIR ${CRMCROOT}/Build)
   message(STATUS "Found CRMC in ${CRMC_LIBRARIES}")
   message(STATUS "CRMC_INCLUDE_DIR = ${CRMC_INCLUDE_DIR}")
   message(STATUS "CRMC_LIBRARY_DIR = ${CRMC_LIBRARY_DIR}")
   message(STATUS "CRMC_BUILD_DIR = ${CRMC_BUILD_DIR}")
else()
   message(STATUS "CRMC library not found; try to set a CRMCROOT environment variable 
           to the base installation path or add -DCRMCROOT = to the cmake command")      
endif()

