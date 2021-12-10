# - build library

if(GEANT4_USE_OPENGL_X11)
  add_definitions(-DG4VIS_USE_OPENGLX)
endif()

if(GEANT4_USE_RAYTRACER_X11)
  add_definitions(-DG4VIS_USE_RAYTRACERX)
endif()

if(GEANT4_USE_XM)
  add_definitions(-DG4VIS_USE_OPENGLXM)
endif()

# library
geant4_add_pymodule(pyG4visualization
  pyG4ASCIITree.cc
  pyG4DAWNFILE.cc
  pyG4HepRepFile.cc
  pyG4OpenGLImmediateX.cc
  pyG4OpenGLImmediateXm.cc
  pyG4OpenGLStoredX.cc
  pyG4OpenGLStoredXm.cc
  pyG4RayTracer.cc
  pyG4RayTracerX.cc
  pyG4VGraphicsSystem.cc
  pyG4VRML2File.cc
  pyG4VisManager.cc
  pymodG4visualization.cc
)

target_link_libraries(pyG4visualization PRIVATE
  G4FR
  G4Tree
  G4RayTracer
  G4VRML
  G4visHepRep
  G4vis_management
)

if(GEANT4_USE_OPENGL_X11 OR GEANT4_USE_XM)
  target_link_libraries(pyG4visualization PRIVATE G4OpenGL)
endif()

install(TARGETS pyG4visualization DESTINATION "${CMAKE_INSTALL_PYTHONDIR}/Geant4")
