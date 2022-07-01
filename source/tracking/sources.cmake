# - G4tracking module build definition

# Define the Geant4 Module.
geant4_add_module(G4tracking
  PUBLIC_HEADERS
    G4AdjointCrossSurfChecker.hh
    G4AdjointSteppingAction.hh
    G4AdjointTrackingAction.hh
    G4RichTrajectory.hh
    G4RichTrajectoryPoint.hh
    G4SmoothTrajectory.hh
    G4SmoothTrajectoryPoint.hh
    G4SteppingManager.hh
    G4SteppingVerbose.hh
    G4SteppingVerboseWithUnits.hh
    G4TrackingManager.hh
    G4TrackingMessenger.hh
    G4Trajectory.hh
    G4TrajectoryPoint.hh
    G4UserSteppingAction.hh
    G4MultiSteppingAction.hh
    G4UserTrackingAction.hh
    G4MultiTrackingAction.hh
    G4VSteppingVerbose.hh
    G4VTrackingManager.hh
    G4VTrajectory.hh
    G4VTrajectoryPoint.hh
    trkgdefs.hh
  SOURCES
    G4AdjointCrossSurfChecker.cc
    G4AdjointSteppingAction.cc
    G4AdjointTrackingAction.cc
    G4RichTrajectory.cc
    G4RichTrajectoryPoint.cc
    G4SmoothTrajectory.cc
    G4SmoothTrajectoryPoint.cc
    G4SteppingManager.cc
    G4SteppingVerbose.cc
    G4SteppingVerboseWithUnits.cc
    G4TrackingManager.cc
    G4TrackingMessenger.cc
    G4Trajectory.cc
    G4TrajectoryPoint.cc
    G4UserSteppingAction.cc
    G4MultiSteppingAction.cc
    G4UserTrackingAction.cc
    G4MultiTrackingAction.cc
    G4VSteppingVerbose.cc
    G4VTrajectory.cc
    G4VTrajectoryPoint.cc
    G4MultiSteppingAction.cc)

geant4_module_compile_definitions(G4tracking PRIVATE G4TRACKING_ALLOC_EXPORT)

geant4_module_link_libraries(G4tracking
  PUBLIC
    G4geometrymng
    G4globman
    G4heprandom
    G4intercoms
    G4navigation
    G4partman
    G4procman
    G4track
    G4volumes
  PRIVATE
    G4detector
    G4graphics_reps)
