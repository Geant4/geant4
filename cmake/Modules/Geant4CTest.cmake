# - Basic setup for testing Geant4 using CTest
# Still rather rough.
# 

#---Deduce the build name--------------------------------------------------------
set(BUILDNAME ${GEANT4_SYSTEM}-${GEANT4_COMPILER}-${CMAKE_BUILD_TYPE})
enable_testing()
include(CTest)
