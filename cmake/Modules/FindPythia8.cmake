# - Locate pythia8 library
# Defines:
#
#  Pythia8_FOUND
#  PYTHIA8_VERSION
#  PYTHIA8_INCLUDE_DIR
#  PYTHIA8_XMLDOC_DIR
#  PYTHIA8_INCLUDE_DIRS (not cached)
#  PYTHIA8_LIBRARY
#  PYTHIA8_hepmcinterface_LIBRARY
#  PYTHIA8_lhapdfdummy_LIBRARY
#  PYTHIA8_LIBRARIES (not cached) : includes 3 libraries above; not to be used if lhapdf is used
set(TEST_PYTHIA8_ROOT_DIR  "" ${PYTHIA8_ROOT_DIR})
IF(TEST_PYTHIA8_ROOT_DIR STREQUAL "")
  IF(DEFINED ENV{PYTHIA8_ROOT_DIR})
    set(PYTHIA8_ROOT_DIR  $ENV{PYTHIA8_ROOT_DIR})
  else()
    set(PYTHIA8_ROOT_DIR  "/usr")
  endif()
endif()

find_path(PYTHIA8_INCLUDE_DIR Pythia.h Pythia8/Pythia.h HINTS  ${PYTHIA8_ROOT_DIR}/include)

find_path(PYTHIA8_XMLDOC_DIR Version.xml HINTS  ${PYTHIA8_ROOT_DIR}/xmldoc  ${PYTHIA8_ROOT_DIR}/share/Pythia8/xmldoc ${PYTHIA8_ROOT_DIR}/share/pythia8-data/xmldoc  ${PYTHIA8_ROOT_DIR}/share/doc/packages/pythia/xmldoc )

if(PYTHIA8_INCLUDE_DIR AND PYTHIA8_XMLDOC_DIR)
  file(READ ${PYTHIA8_XMLDOC_DIR}/Version.xml versionstr)
  string(REGEX REPLACE ".*Pythia:versionNumber.*default.*[0-9][.]([0-9]+).*" "\\1" PYTHIA8_VERSION "${versionstr}")
  set(PYTHIA8_VERSION "8.${PYTHIA8_VERSION}")
  find_library(PYTHIA8_LIBRARY NAMES pythia8 Pythia8 HINTS ${PYTHIA8_ROOT_DIR}/lib ${PYTHIA8_ROOT_DIR}/lib64)
  find_library(PYTHIA8_lhapdfdummy_LIBRARY NAMES lhapdfdummy  HINTS ${PYTHIA8_ROOT_DIR}/lib ${PYTHIA8_ROOT_DIR}/lib64)
  set(PYTHIA8_INCLUDE_DIRS ${PYTHIA8_INCLUDE_DIR} ${PYTHIA8_INCLUDE_DIR}/Pythia8 ${PYTHIA8_INCLUDE_DIR}/Pythia8Plugins )
  set(PYTHIA8_LIBRARIES ${PYTHIA8_LIBRARY})
  if(PYTHIA8_VERSION VERSION_LESS 8.200)
    #Is this library needed?
    set(PYTHIA8_LIBRARIES ${PYTHIA8_LIBRARY} ${PYTHIA8_lhapdfdummy_LIBRARY})
  endif()
  find_file(resHEPMC3 HepMC3.h PATHS  ${PYTHIA8_INCLUDE_DIRS} NO_DEFAULT_PATH)
  if (resHEPMC3)
    set(Pythia8_HEPMC3_FOUND TRUE)
  endif()
endif()

# handle the QUIETLY and REQUIRED arguments and set Pythia8_FOUND to TRUE if
# all listed variables are TRUE

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Pythia8 REQUIRED_VARS PYTHIA8_INCLUDE_DIR PYTHIA8_LIBRARIES PYTHIA8_XMLDOC_DIR VERSION_VAR PYTHIA8_VERSION HANDLE_COMPONENTS)

mark_as_advanced(Pythia8_FOUND PYTHIA8_INCLUDE_DIR PYTHIA8_LIBRARY PYTHIA8_LIBRARIES PYTHIA8_XMLDOC_DIR)

