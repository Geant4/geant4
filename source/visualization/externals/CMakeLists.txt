#------------------------------------------------------------------------------
# Module : 
# Package: Geant4.src.G4visualization.
#------------------------------------------------------------------------------

#
# GL2PS built for GL and Inventor drivers
#
if(GEANT4_USE_OPENGL OR GEANT4_USE_INVENTOR)
  add_subdirectory(gl2ps)
endif()
