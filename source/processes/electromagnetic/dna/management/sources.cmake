#------------------------------------------------------------------------------
# sources.cmake
# Module : G4emlowenergy
# Package: Geant4.src.G4processes.G4electromagnetic.G4emlowenergy
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 101354 2016-11-15 08:27:51Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/magneticfield/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/cuts/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/utils/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4emdna-man
    HEADERS
        AddClone_def.hh
        G4AllITFinder.hh
        G4CTCounter.hh
        G4FastList.hh
        G4FastList.icc
        G4IosFlagsSaver.hh
        G4ITBox.hh
        G4IT.hh
        G4ITGun.hh
        G4ITFinder.hh
        G4ITFinder.icc
        G4ITLeadingTracks.hh
        G4ITModelHandler.hh
        G4ITModelManager.hh
        G4ITModelProcessor.hh
        G4ITMultiNavigator.hh
        G4ITNavigator.hh
        G4ITNavigator1.hh
        G4ITNavigator1.icc
        G4ITNavigator2.hh
        G4ITNavigator2.icc
        G4ITPathFinder.hh
        G4ITReaction.hh
        G4ITReactionChange.hh
        G4ITReactionTable.hh
        G4ITSafetyHelper.hh
        G4ITStepProcessor.hh
        G4ITStepStatus.hh
        G4ITSteppingVerbose.hh
        G4ITTrackHolder.hh
        G4ITTrackingInteractivity.hh
        G4ITTrackingManager.hh
        G4ITTransportation.hh
        G4ITTransportation.icc
        G4ITTransportationManager.hh
        G4ITTransportationManager.icc
        G4ITType.hh
        G4ITStepProcessorState_Lock.hh
        G4KDMap.hh
        G4KDNode.hh
        G4KDNode.icc
        G4KDTree.hh
        G4KDTree.icc
        G4KDTreeResult.hh
        G4ManyFastLists.hh
        G4ManyFastLists.icc
        G4MemStat.hh
        G4ReferenceCast.hh
        G4memory.hh
        G4Scheduler.hh
        G4SchedulerMessenger.hh
        G4TrackingInformation.hh
        G4TrackList.hh
        G4TrackState.hh
        G4UserTimeStepAction.hh
        G4VITDiscreteProcess.hh
        G4VITProcess.hh
        G4VITReactionProcess.hh
        G4VITRestDiscreteProcess.hh
        G4VITRestProcess.hh
        G4VITStepModel.hh
        G4VITSteppingVerbose.hh
        G4VITTimeStepComputer.hh
        G4VITTrackHolder.hh
        G4VScheduler.hh
    SOURCES
        G4AllITFinder.cc
        G4ITBox.cc
        G4IT.cc
        G4ITGun.cc
        G4ITFinder.cc
        G4ITLeadingTracks.cc
        G4ITModelHandler.cc
        G4ITModelManager.cc
        G4ITModelProcessor.cc
        G4ITMultiNavigator.cc
        G4ITNavigator1.cc
        G4ITNavigator2.cc
        G4ITNavigatorState2.cc
        G4ITPathFinder.cc
        G4ITReaction.cc
        G4ITReactionChange.cc
        G4ITReactionTable.cc
        G4ITSafetyHelper.cc
        G4ITStepProcessor2.cc
        G4ITStepProcessor.cc
        G4ITSteppingVerbose.cc
        G4ITTrackHolder.cc
        G4ITTrackingInteractivity.cc
        G4ITTrackingManager.cc
        G4ITTransportation.cc
        G4ITTransportationManager.cc
        G4ITType.cc
        G4KDMap.cc
        G4KDNode.cc
        G4KDTree.cc
        G4KDTreeResult.cc
        G4MemStat.cc
        G4Scheduler.cc
        G4SchedulerMessenger.cc
        G4TrackingInformation.cc
        G4TrackList.cc
        G4TrackState.cc
        G4UserTimeStepAction.cc
        G4VITDiscreteProcess.cc
        G4VITProcess.cc
        G4VITReactionProcess.cc
        G4VITRestDiscreteProcess.cc
        G4VITRestProcess.cc
        G4VITStepModel.cc
        G4VITSteppingVerbose.cc
        G4VITTimeStepComputer.cc
        G4VITTrackHolder.cc
        G4VScheduler.cc
    GRANULAR_DEPENDENCIES
        G4detector
        G4geometrymng
        G4magneticfield
        G4globman
        G4intercoms
        G4materials
        G4navigation
        G4volumes
        G4partman
        G4procman
        G4cuts
        G4emutils
    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4intercoms
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

