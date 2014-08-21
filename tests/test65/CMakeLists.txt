cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test65)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

#-----------------------------------------------------------------------
# Download the required dataset and update testing environment to use it
#
include(ExternalProject)

#temporary replace by local download, original file has symlink, local not.
#URL ftp://gdo-nuclear.ucllnl.org/pub/G4LEND/ENDF.B-VII.0.tar.gz

ExternalProject_Add(test65-data
  PREFIX Data/tmp
  SOURCE_DIR Data/g4lend
  URL http://geant4.cern.ch/support/source/ENDF.B-VII.0.tar.gz
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  )

#temporary replaced by local download:
#URL ftp://gdo-nuclear.ucllnl.org/pub/G4LEND/ENDF.B-VII.0.tar.gz
list(APPEND GEANT4_TEST_ENVIRONMENT "G4LENDDATA=${PROJECT_BINARY_DIR}/Data/g4lend")

GEANT4_EXECUTABLE(test65 test65.cc src/*.cc)
GEANT4_EXECUTABLE(test65_hadronic_exerciser test65_hadronic_exerciser.cc)

# Ensure data is built first
add_dependencies(test65 test65-data)

add_dependencies(test65 test65_hadronic_exerciser) # Artificial to force serialization (problem on windows)

#-----------------------------------------------------------------------
# Test definitions
#
GEANT4_ADD_TEST(test65_hadronic_exerciser COMMAND test65_hadronic_exerciser
                                          OUTPUT  test65.in
                                          BUILD test65_hadronic_exerciser)
GEANT4_ADD_TEST(test65 COMMAND test65 test65.in
                       BUILD test65
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                       DEPENDS test65_hadronic_exerciser)
