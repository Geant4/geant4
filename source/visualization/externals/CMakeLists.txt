# - gl2ps only needed if OGL drivers enabled
if(GEANT4_USE_OPENGL OR GEANT4_USE_INVENTOR)
  add_subdirectory(gl2ps)
endif()
