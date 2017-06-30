#------------------------------------------------------------------------------
# sources.cmake
# Module : G4globman
# Package: Geant4.src.G4global.G4globman
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 103592 2017-04-19 08:09:03Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4globman
    HEADERS
        globals.hh
        templates.hh
        tls.hh
        windefs.hh
        G4Allocator.hh
        G4strstreambuf.icc
        G4AllocatorPool.hh
        G4AllocatorList.hh
        G4ApplicationState.hh
        G4AutoLock.hh
        G4DataVector.hh
        G4DataVector.icc
        G4ErrorPropagatorData.hh
        G4ErrorPropagatorData.icc
        G4Evaluator.hh
        G4ExceptionSeverity.hh
        G4Exp.hh
        G4FPEDetection.hh
        G4FastVector.hh
        G4GeometryTolerance.hh
        G4Log.hh
        G4LPhysicsFreeVector.hh
        G4OrderedTable.hh
        G4PhysicalConstants.hh
        G4PhysicsFreeVector.hh
        G4PhysicsLinearVector.hh
        G4PhysicsLnVector.hh
        G4PhysicsLogVector.hh
        G4PhysicsModelCatalog.hh
        G4PhysicsOrderedFreeVector.hh
        G4PhysicsTable.hh
        G4PhysicsTable.icc
        G4PhysicsVector.hh
        G4PhysicsVector.icc
        G4PhysicsVectorType.hh
        G4Physics2DVector.hh
        G4Physics2DVector.icc
        G4Pow.hh
        G4ReferenceCountedHandle.hh
        G4RotationMatrix.hh
        G4SIunits.hh
        G4SliceTimer.hh
        G4SliceTimer.icc
        G4StateManager.hh
        G4StateManager.icc
        G4String.hh
        G4String.icc
        G4SystemOfUnits.hh
        G4Threading.hh
        G4ThreeVector.hh
        G4Timer.hh
        G4Timer.icc
        G4Tokenizer.hh
        G4TwoVector.hh
        G4Types.hh
        G4UnitsTable.hh
        G4UnitsTable.icc
        G4UserLimits.hh
        G4UserLimits.icc
        G4VExceptionHandler.hh
        G4VNotifier.hh
        G4VStateDependent.hh
        G4Version.hh
        G4coutDestination.hh
        G4ios.hh
        G4strstreambuf.hh
        G4MTcoutDestination.hh
        G4CacheDetails.hh
        G4Cache.hh
        G4ThreadLocalSingleton.hh
        G4AutoDelete.hh
        G4TWorkspacePool.hh
        G4MTBarrier.hh
        G4coutFormatters.hh
        G4MulticoutDestination.hh
        G4LockcoutDestination.hh
        G4MasterForwardcoutDestination.hh
        G4FilecoutDestination.hh
        G4BuffercoutDestination.hh
    SOURCES
        G4Allocator.cc
        G4AllocatorPool.cc
        G4AllocatorList.cc
        G4DataVector.cc
        G4ErrorPropagatorData.cc
        G4Exception.cc
        G4GeometryTolerance.cc
        G4LPhysicsFreeVector.cc
        G4OrderedTable.cc
        G4PhysicsFreeVector.cc
        G4PhysicsLinearVector.cc
        G4PhysicsLogVector.cc
        G4PhysicsModelCatalog.cc
        G4PhysicsOrderedFreeVector.cc
        G4PhysicsTable.cc
        G4PhysicsVector.cc
        G4Physics2DVector.cc
        G4Pow.cc
        G4ReferenceCountedHandle.cc
        G4SliceTimer.cc
        G4StateManager.cc
        G4Threading.cc
        G4Timer.cc
        G4UnitsTable.cc
        G4VExceptionHandler.cc
        G4VNotifier.cc
        G4VStateDependent.cc
        G4coutDestination.cc
        G4ios.cc
        G4MTcoutDestination.cc
        G4CacheDetails.cc
        G4MTBarrier.cc
        G4coutFormatters.cc
        G4LockcoutDestination.cc
        G4MasterForwardcoutDestination.cc
        G4FilecoutDestination.cc
        G4BuffercoutDestination.cc
    GRANULAR_DEPENDENCIES
    GLOBAL_DEPENDENCIES
    LINK_LIBRARIES
        ${CLHEP_LIBRARIES}
)

# List any source specific properties here

