#------------------------------------------------------------------------------
# Module : G4RayTracer
# Package: Geant4.src.G4visualization.G4RayTracer
#------------------------------------------------------------------------------

#
# Module has optional sources
#
# List those always built
set(G4VIS_RAYTRACER_MODULE_HEADERS
  G4RTJpeg.hh
  G4RTJpegCoder.hh
  G4RTJpegCoderTables.hh
  G4RTJpegMaker.hh
  G4RTMessenger.hh
  G4RTOutBitStream.hh
  G4RTPrimaryGeneratorAction.hh
  G4RTRun.hh
  G4RTRunAction.hh
  G4RTSimpleScanner.hh
  G4RTSteppingAction.hh
  G4RTTrackingAction.hh
  G4RTWorkerInitialization.hh
  G4RayTracer.hh
  G4RayTracerFeatures.hh
  G4RayTracerSceneHandler.hh
  G4RayTracerViewer.hh
  G4RayTrajectory.hh
  G4RayTrajectoryPoint.hh
  G4TheMTRayTracer.hh
  G4TheRayTracer.hh
  G4VFigureFileMaker.hh
  G4VRTScanner.hh)

set(G4VIS_RAYTRACER_MODULE_SOURCES
  G4RTJpegCoder.cc
  G4RTJpegMaker.cc
  G4RTMessenger.cc
  G4RTOutBitStream.cc
  G4RTPrimaryGeneratorAction.cc
  G4RTRun.cc
  G4RTRunAction.cc
  G4RTSimpleScanner.cc
  G4RTSteppingAction.cc
  G4RTTrackingAction.cc
  G4RTWorkerInitialization.cc
  G4RayTracer.cc
  G4RayTracerSceneHandler.cc
  G4RayTracerViewer.cc
  G4RayTrajectory.cc
  G4RayTrajectoryPoint.cc
  G4TheMTRayTracer.cc
  G4TheRayTracer.cc
  G4VRTScanner.cc)

set(G4VIS_RAYTRACER_MODULE_LINK_LIBRARIES )

#
# X11 RayTracer only if selected
#
if(GEANT4_USE_RAYTRACER_X11)
  list(APPEND G4VIS_RAYTRACER_MODULE_HEADERS
    G4RayTracerX.hh
    G4RayTracerXViewer.hh
    G4RTXScanner.hh)

  list(APPEND G4VIS_RAYTRACER_MODULE_SOURCES
    G4RayTracerX.cc
    G4RayTracerXViewer.cc
    G4RTXScanner.cc)

  #
  # Add source properties and additional LINK_LIBRARIES here
  #
  # Must use G4VIS_BUILD_RAYTRACERX_DRIVER define
  add_definitions(-DG4VIS_BUILD_RAYTRACERX_DRIVER)

  # The X11 Libraries
  list(APPEND G4VIS_RAYTRACER_MODULE_LINK_LIBRARIES X11::SM X11::ICE X11::X11 X11::Xext X11::Xmu)
endif()

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4RayTracer
  HEADERS
    ${G4VIS_RAYTRACER_MODULE_HEADERS}
  SOURCES
    ${G4VIS_RAYTRACER_MODULE_SOURCES}
  GRANULAR_DEPENDENCIES
    G4bosons
    G4cuts
    G4detector
    G4digits
    G4event
    G4geometrymng
    G4globman
    G4graphics_reps
    G4hits
    G4intercoms
    G4materials
    G4modeling
    G4navigation
    G4partman
    G4procman
    G4track
    G4tracking
    G4vis_management
    G4volumes
  GLOBAL_DEPENDENCIES
    G4digits_hits
    G4event
    G4geometry
    G4global
    G4graphics_reps
    G4intercoms
    G4materials
    G4modeling
    G4particles
    G4processes
    G4track
    G4tracking
    G4vis_management
  LINK_LIBRARIES
    ${G4VIS_RAYTRACER_MODULE_LINK_LIBRARIES}
)

# List any source specific properties here
