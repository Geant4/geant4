# - G4RayTracer module build definition

# Module has optional sources
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

# X11 RayTracer only if selected
if(GEANT4_USE_RAYTRACER_X11)
  list(APPEND G4VIS_RAYTRACER_MODULE_HEADERS
    G4RayTracerX.hh
    G4RayTracerXViewer.hh
    G4RTXScanner.hh)

  list(APPEND G4VIS_RAYTRACER_MODULE_SOURCES
    G4RayTracerX.cc
    G4RayTracerXViewer.cc
    G4RTXScanner.cc)

  # The X11 Libraries
  list(APPEND G4VIS_RAYTRACER_MODULE_LINK_LIBRARIES X11::SM X11::ICE X11::X11 X11::Xext X11::Xmu)
endif()

# Define the Geant4 Module.
geant4_add_module(G4RayTracer
  PUBLIC_HEADERS ${G4VIS_RAYTRACER_MODULE_HEADERS}
  SOURCES ${G4VIS_RAYTRACER_MODULE_SOURCES})

geant4_module_link_libraries(G4RayTracer
  PUBLIC
    G4globman
    G4graphics_reps
    G4modeling
    G4run
    G4hits
    G4track
    G4intercoms
    G4tracking
    G4vis_management
    G4hepgeometry
    ${G4VIS_RAYTRACER_MODULE_LINK_LIBRARIES}
  PRIVATE
    G4event
    G4bosons
    G4scoring
    G4partman
    G4procman
    G4cuts
    G4tasking
    G4detector
    G4navigation)
