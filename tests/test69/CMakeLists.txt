cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test69)
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

# write the location of the macros in the Tst69MacroLocation.hh header
configure_file(${PROJECT_SOURCE_DIR}/include/Tst69MacroLocation.hh.in
               ${PROJECT_BINARY_DIR}/include/Tst69MacroLocation.hh)
include_directories(BEFORE ${PROJECT_BINARY_DIR}/include)

GEANT4_EXECUTABLE(test69 test69.cc src/*.cc)

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

set(TEST69_AUXILIARY_MACROS
  vis.mac
  )

foreach(_macro ${TEST69_MACROS} ${TEST69_SHORT_MACROS} ${TEST69_AUXILIARY_MACROS})
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
GEANT4_ADD_TEST(test69-build BUILD test69 LABELS Nightly Continuous) 

foreach(_macro ${TEST69_MACROS})
  get_filename_component(_macro_name ${_macro} NAME_WE)

  foreach(_list ${TEST69_PHYSLISTS})

    GEANT4_ADD_TEST(test69-${_macro_name}-${_list}
                    COMMAND test69 ${PROJECT_BINARY_DIR}/macros/${_macro}
                    DEPENDS test69-build
                    LABELS Nightly Continuous
                    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} "TEST69_PHYSLIST=${_list}")

  endforeach()
endforeach()

foreach(_macro ${TEST69_SHORT_MACROS})
  get_filename_component(_macro_name ${_macro} NAME_WE)

  foreach(_list ${TEST69_HP_PHYSLISTS})

    GEANT4_ADD_TEST(test69-${_macro_name}-${_list}
                    COMMAND test69 ${PROJECT_BINARY_DIR}/macros/${_macro}
                    DEPENDS test69-build
                    LABELS Nightly Continuous
                    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} "TEST69_PHYSLIST=${_list}")

  endforeach()
endforeach()

GEANT4_ADD_TEST(test69-p_40GeV_mix-FTFP_INCLXX_ABLA
  COMMAND test69 ${PROJECT_BINARY_DIR}/macros/p_40GeV_mix.in
  DEPENDS test69-build
  LABELS Nightly Continuous
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} "TEST69_PHYSLIST=FTFP_INCLXX" "TEST69_USE_ABLA=1")

