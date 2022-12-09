# - Builds for always-on vis libraries
add_subdirectory(FukuiRenderer)
add_subdirectory(HepRep)
add_subdirectory(RayTracer)
add_subdirectory(Tree)
add_subdirectory(VRML)
add_subdirectory(gMocren)
add_subdirectory(management)
add_subdirectory(modeling)
add_subdirectory(ToolsSG)

# OpenGL is optional depending on user selection
if(GEANT4_USE_OPENGL)
  add_subdirectory(OpenGL)
endif()

# OpenInventor is optional depending on user selection
if(GEANT4_USE_INVENTOR)
  add_subdirectory(OpenInventor)
endif()

# Qt3D is available if Qt is used and > 5.15
if(GEANT4_USE_QT3D)
  add_subdirectory(Qt3D)
endif()

# VTK is optional depending on user selection
if(GEANT4_USE_VTK)
  add_subdirectory(Vtk)
endif()