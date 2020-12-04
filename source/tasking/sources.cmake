#------------------------------------------------------------------------------
# sources.cmake
# Module : G4run
# Package: Geant4.src.G4run
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 110666 2018-06-06 15:09:06Z jmadsen $
#
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
GEANT4_DEFINE_MODULE(NAME G4tasking
    HEADERS
        G4RunManagerFactory.hh
        G4TaskRunManager.hh
        G4TaskRunManagerKernel.hh
        G4UserTaskInitialization.hh
        G4UserTaskThreadInitialization.hh
        G4Task.hh
        G4TaskGroup.hh
        G4TBBTaskGroup.hh
        G4TaskSingletonDelegator.hh
        G4ThreadData.hh
        G4ThreadPool.hh
        G4VTask.hh
        G4VTaskGroup.hh
        G4TaskManager.hh
        G4UserTaskQueue.hh
        G4VUserTaskQueue.hh
        G4WorkerTaskRunManager.hh
        G4WorkerTaskRunManagerKernel.hh
        taskdefs.hh
   SOURCES
        G4RunManagerFactory.cc
        G4TaskRunManager.cc
        G4TaskRunManagerKernel.cc
        G4UserTaskInitialization.cc
        G4UserTaskThreadInitialization.cc
        G4WorkerTaskRunManager.cc
        G4WorkerTaskRunManagerKernel.cc
    GRANULAR_DEPENDENCIES
        G4cuts
        G4decay
        G4detector
        G4detutils
        G4digits
        G4emutils
        G4event
        G4geombias
        G4geometrymng
        G4globman
        G4graphics_reps
        G4hepnumerics
        G4hits
        G4intercoms
        G4magneticfield
        G4materials
        G4navigation
        G4partman
        G4procman
        G4run
        G4scoring
        G4track
        G4tracking
        G4transportation
        G4volumes
        G4specsolids
    GLOBAL_DEPENDENCIES
        G4digits_hits
        G4event
        G4geometry
        G4global
        G4graphics_reps
        G4intercoms
        G4materials
        G4particles
        G4processes
        G4run
        G4track
        G4tracking
    LINK_LIBRARIES
        ${timemory_LIBRARIES}
        ${TBB_LIBRARIES}
        ${PTL_LIBRARIES}
)

# List any source specific properties here

