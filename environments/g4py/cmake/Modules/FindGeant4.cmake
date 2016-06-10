# - Find Geant4 library
# This module sets up Geant4 information
# It defines:
# GEANT4_FOUND               If the Geant4 is found
# GEANT4_INCLUDE_DIR         PATH to the include directory
# GEANT4_LIBRARY_DIR         PATH to the library directory
# GEANT4_LIBRARIES           Most common libraries
# GEANT4_LIBRARIES_WITH_VIS  Most common libraries with visualization

find_program(GEANT4_CONFIG NAMES geant4-config
             PATHS $ENV{GEANT4_INSTALL}/bin
                   ${GEANT4_INSTALL}/bin
                   /usr/local/bin /opt/local/bin)

if(GEANT4_CONFIG)
  set(GEANT4_FOUND TRUE)

  execute_process(COMMAND ${GEANT4_CONFIG} --prefix
                  OUTPUT_VARIABLE GEANT4_PREFIX
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${GEANT4_CONFIG} --version
                  OUTPUT_VARIABLE GEANT4_VERSION
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(COMMAND ${GEANT4_CONFIG} --has-feature gdml
                  OUTPUT_VARIABLE _TMP)
  if (_TMP MATCHES "yes")
    set(GEANT4_HAS_GDML TRUE)
  else()
    set(GEANT4_HAS_GDML FALSE)
  endif()

  execute_process(COMMAND ${GEANT4_CONFIG} --has-feature opengl-x11
                  OUTPUT_VARIABLE _TMP)
  if (_TMP MATCHES "yes")
    set(GEANT4_HAS_OPENGL TRUE)
  else()
    set(GEANT4_HAS_OPENGL FALSE)
  endif()

  execute_process(COMMAND ${GEANT4_CONFIG} --has-feature raytracer-x11
                  OUTPUT_VARIABLE _TMP)
  if (_TMP MATCHES "yes")
    set(GEANT4_HAS_RAYTRACER_X11 TRUE)
  else()
    set(GEANT4_HAS_RAYTRACER_X11 FALSE)
  endif()

  execute_process(COMMAND ${GEANT4_CONFIG} --has-feature motif
                  OUTPUT_VARIABLE _TMP)
  if (_TMP MATCHES "yes")
    set(GEANT4_HAS_MOTIF TRUE)
  else()
    set(GEANT4_HAS_MOTIF FALSE)
  endif()

  message(STATUS "Found Geant4: ${GEANT4_PREFIX} (${GEANT4_VERSION})")

else()
  set(GEANT4_FOUND FALSE)
  message(SEND_ERROR "NOT Found Geant4: set GEANT4_INSTALL env.")

endif()

set(GEANT4_INCLUDE_DIR ${GEANT4_PREFIX}/include/Geant4)
set(GEANT4_LIBRARY_DIR ${GEANT4_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(GEANT4_LIBRARIES  G4interfaces G4persistency G4analysis
                      G4error_propagation G4readout G4physicslists
                      G4run G4event G4tracking G4parmodels G4processes
                      G4digits_hits G4track G4particles G4geometry
                      G4materials G4graphics_reps G4intercoms
                      G4global G4clhep)

set(GEANT4_LIBRARIES_WITH_VIS
                      G4OpenGL G4gl2ps G4Tree G4FR G4GMocren G4visHepRep
                      G4RayTracer G4VRML G4vis_management G4modeling
                      G4interfaces G4persistency G4analysis
                      G4error_propagation G4readout G4physicslists
                      G4run G4event G4tracking G4parmodels G4processes
                      G4digits_hits G4track G4particles G4geometry
                      G4materials G4graphics_reps G4intercoms
                      G4global G4clhep)
