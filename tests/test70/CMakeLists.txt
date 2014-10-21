cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test70)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test70 test70.cc src/*.cc)


# ----------------------------------------------------------------------------------------------------
# --- Copy MAC files 
# ----------------------------------------------------------------------------------------------------

#configure_file(short.mac short.mac COPYONLY) # test expects file to be local
#configure_file(short-mt.mac short-mt.mac COPYONLY) # test expects file to be local
#configure_file(long.mac long.mac COPYONLY) # test expects file to be local

#-------------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
        
file(GLOB MAC_FILES ${PROJECT_SOURCE_DIR}/*.mac ${PROJECT_SOURCE_DIR}/*.in)

foreach(_script ${MAC_FILES})
  configure_file(
    ${_script}
    ${PROJECT_BINARY_DIR}/.
    COPYONLY
    )
endforeach()


# ----------------------------------------------------------------------------------------------------
# --- Test definitions
# ----------------------------------------------------------------------------------------------------

GEANT4_ADD_TEST(test70-build BUILD test70
                       LABELS Nightly Continuous
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})


GEANT4_ADD_TEST(test70 COMMAND test70 -mac ${CMAKE_CURRENT_SOURCE_DIR}/short.mac
                       DEPENDS test70-build
                       LABELS Nightly Continuous
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})


GEANT4_ADD_TEST(test70-mt COMMAND test70 -mt 2 -mac ${CMAKE_CURRENT_SOURCE_DIR}/short-mt.mac
                       DEPENDS test70-build
                       LABELS Nightly Continuous
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

# ----------------------------------------------------------------------------------------------------
# --- Nightly ONLY
# ----------------------------------------------------------------------------------------------------

#GEANT4_ADD_TEST(test70-long COMMAND test70 -mac ${CMAKE_CURRENT_SOURCE_DIR}/long.mac
#                       DEPENDS test70-build
#                       LABELS Nightly
#                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})


#GEANT4_ADD_TEST(test70-long-mt COMMAND test70 -mt 10 -mac ${CMAKE_CURRENT_SOURCE_DIR}/long.mac
#                       DEPENDS test70-build
#                       LABELS Nightly
#                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

