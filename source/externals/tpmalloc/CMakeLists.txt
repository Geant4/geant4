#-----------------------------------------------------------------------
# CMakeLists.txt
# Module : G4tpmalloc
# Package: Geant4.src.G4externals.G4tpmalloc
#
# CMakeLists.txt for G4tpmalloc. We do not use the usual Geant4 sources.cmake
# approach because G4tpmalloc requires a custom install solution for its
# headers. These are in a nested structure so this needs to be replicated
# when we install, i.e:
# <INCLUDE>
# +- geant4
#    +- tpmalloc
#       +- mymalloc.h
#       ...
#
# Created on: 17/04/2013
#
# $Date$
# $Revision$
# $Author$
#
#------------------------------------------------------------------------------

#-----------------------------------------------------------------------
# Force the compiler flags to always use -O3
# Handle both build and configuration types
set(CMAKE_C_FLAGS "-O3")

string(TOUPPER "${CMAKE_BUILD_TYPE}" TPMALLOC_BUILD_TYPE)
set(CMAKE_C_FLAGS_${TPMALLOC_BUILD_TYPE} "")

foreach(_type ${CMAKE_CONFIGURATION_TYPES})
  string(TOUPPER "${_type}" _upper_type)
  set(CMAKE_C_FLAGS_${TPMALLOC_BUILD_TYPE} "")
endforeach()

#----------------------------------------------------------------------------
# Includes for this modules
#
include_directories(include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)

#----------------------------------------------------------------------------
# List the headers and sources
#
# - Headers
set(G4tpmalloc_HEADERS
    include/tpmalloc/mymalloc.h
    include/tpmalloc/tpmallocstub.h
)

# - Sources
set(G4tpmalloc_SOURCES
    src/hjmalloc.cc
    src/mymalloc.cc
)

#----------------------------------------------------------------------------
# Now add the library targets
#
include(Geant4MacroLibraryTargets)
GEANT4_LIBRARY_TARGET(NAME G4tpmalloc SOURCES ${G4tpmalloc_SOURCES} ${G4tpmalloc_HEADERS})

#----------------------------------------------------------------------------
# Ensure the header directory gets added to the list of ones to export
#
set_property(GLOBAL APPEND PROPERTY GEANT4_BUILDTREE_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

#----------------------------------------------------------------------------
# Install the headers, retaining their nested structure
#
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    COMPONENT Development
    PATTERN ".svn" EXCLUDE
    PATTERN "*.h")
