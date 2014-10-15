cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test69)
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

find_package(ROOT QUIET)
if(NOT ROOT_FOUND)
  message(STATUS "G4 TESTS: ROOT package not found. --> ROOT tally disabled in test69")
else()
  add_definitions(-DTEST69_HAS_ROOT)
  include_directories(${ROOT_INCLUDE_DIR})
endif()

# write the location of the macros in the Tst69MacroLocation.hh header
configure_file(${PROJECT_SOURCE_DIR}/include/Tst69MacroLocation.hh.in
               ${PROJECT_BINARY_DIR}/include/Tst69MacroLocation.hh)
include_directories(BEFORE ${PROJECT_BINARY_DIR}/include)

GEANT4_EXECUTABLE(test69 test69.cc src/*.cc LIBRARIES ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all macros to the build directory, i.e. the directory in which we
# build test69. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(TEST69_MACROS
  p_40GeV_mix.in
  C12_12GeV_mix.in
  )

set(TEST69_SHORT_MACROS
  p_40GeV_mix_short.in
  C12_12GeV_mix_short.in
  )

set(TEST69_LONG_MACROS
  p_40GeV_mix_long.in
  )

set(TEST69_AUXILIARY_MACROS
  vis.mac
  )

foreach(_macro ${TEST69_MACROS} ${TEST69_SHORT_MACROS} ${TEST69_LONG_MACROS} ${TEST69_AUXILIARY_MACROS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/macros/${_macro}
    ${PROJECT_BINARY_DIR}/macros/${_macro}
    COPYONLY
    )
endforeach()

set(TEST69_PHYSLISTS
  QGSP_INCLXX
  FTFP_INCLXX
  )

set(TEST69_HP_PHYSLISTS
  QGSP_INCLXX_HP
  FTFP_INCLXX_HP
  )

#---Test definitions----------------------------------------------------------------------------------
# We define a test69-build test that only builds the executable.
# We also define a normal test for each existing .in file.
GEANT4_ADD_TEST(test69-build BUILD test69 LABELS Nightly Continuous PhysicsChecks) 

foreach(_macro ${TEST69_MACROS})
  get_filename_component(_macro_name ${_macro} NAME_WE)

  foreach(_list ${TEST69_PHYSLISTS})

    GEANT4_ADD_TEST(test69-${_macro_name}-${_list}
                    COMMAND test69 ${PROJECT_BINARY_DIR}/macros/${_macro}
                    DEPENDS test69-build
                    LABELS Nightly Continuous
                    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} "TEST69_PHYSLIST=${_list}"
                    TIMEOUT 3600)

  endforeach()
endforeach()

foreach(_macro ${TEST69_SHORT_MACROS})
  get_filename_component(_macro_name ${_macro} NAME_WE)

  foreach(_list ${TEST69_HP_PHYSLISTS})

    GEANT4_ADD_TEST(test69-${_macro_name}-${_list}
                    COMMAND test69 ${PROJECT_BINARY_DIR}/macros/${_macro}
                    DEPENDS test69-build
                    LABELS Nightly Continuous
                    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} "TEST69_PHYSLIST=${_list}"
                    TIMEOUT 3600)

  endforeach()
endforeach()


# PhysicsChecks tests

GEANT4_ADD_TEST(test69-p_40GeV_mix-FTFP_INCLXX-physics
  COMMAND test69 ${PROJECT_BINARY_DIR}/macros/p_40GeV_mix_long.in
  DEPENDS test69-build
  LABELS PhysicsChecks
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} "TEST69_PHYSLIST=FTFP_INCLXX" "TEST69_TALLY=analysis"
  TIMEOUT 3600)

FIND_PACKAGE(StatTest QUIET) 

if(STATTEST_FOUND)      
  STATTEST_ADD_TEST(test69-p_40GeV_mix-FTFP_INCLXX-physicschecks
    G4TEST test69-p_40GeV_mix-FTFP_INCLXX-physics
    CONFIG "${PROJECT_SOURCE_DIR}/stattest/conf.qa"
    REFERENCE "${GEANT4_TEST_REFERENCES_URL}/tests/test69/test69_tally_FTFP_INCLXX-current.root"
    INPUT "${PROJECT_BINARY_DIR}/test69_tally_FTFP_INCLXX.root"
    IMG "test69_tally_FTFP_INCLXX-physicschecks.pdf"
    LABELS PhysicsChecks
    )
endif(STATTEST_FOUND)

