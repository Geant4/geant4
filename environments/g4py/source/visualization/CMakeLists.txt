# - build library

if(GEANT4_HAS_OPENGL)
  add_definitions(-DG4VIS_USE_OPENGLX)
endif()

if(GEANT4_HAS_RAYTRACER_X11)
  add_definitions(-DG4VIS_USE_RAYTRACERX)
endif()

if(GEANT4_HAS_MOTIF)
  add_definitions(-DG4VIS_USE_OPENGLXM)
endif()

# library
set(_TARGET pyG4visualization)
add_library(
  ${_TARGET} SHARED
  pyG4ASCIITree.cc
  pyG4DAWNFILE.cc
  pyG4HepRep.cc
  pyG4HepRepFile.cc
  pyG4OpenGLImmediateX.cc
  pyG4OpenGLImmediateXm.cc
  pyG4OpenGLStoredX.cc
  pyG4OpenGLStoredXm.cc
  pyG4RayTracer.cc
  pyG4RayTracerX.cc
  pyG4VGraphicsSystem.cc
  pyG4VRML1File.cc
  pyG4VRML2File.cc
  pyG4VisManager.cc
  pymodG4visualization.cc
)

set_target_properties(${_TARGET} PROPERTIES PREFIX "")
set_target_properties(${_TARGET} PROPERTIES OUTPUT_NAME "G4visualization")
set_target_properties(${_TARGET} PROPERTIES SUFFIX ".so")
set_target_properties(${_TARGET}
                      PROPERTIES INSTALL_RPATH
                      ${GEANT4_LIBRARY_DIR}
                      BUILD_WITH_INSTALL_RPATH TRUE)

target_link_libraries (${_TARGET}
                      ${GEANT4_LIBRARIES_WITH_VIS} boost_python
                      ${PYTHON_LIBRARIES})

# install
install(TARGETS ${_TARGET} LIBRARY DESTINATION ${G4MODULES_INSTALL_DIR})
