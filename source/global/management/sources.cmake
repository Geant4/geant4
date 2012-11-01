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
# $Id: sources.cmake,v 1.1 2010-09-29 18:44:49 bmorgan Exp $
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
        G4strstreambuf.hh
        G4Allocator.hh
        G4strstreambuf.icc
        G4AllocatorPool.hh
        G4ApplicationState.hh
        G4DataVector.hh
        G4DataVector.icc
        G4ErrorPropagatorData.hh
        G4ErrorPropagatorData.icc
        G4Evaluator.hh
        G4ExceptionSeverity.hh
        G4FPEDetection.hh
        G4FastVector.hh
        G4GeometryTolerance.hh
        G4LPhysicsFreeVector.hh
        G4LPhysicsFreeVector.icc
        G4OrderedTable.hh
        G4PhysicalConstants.hh
        G4PhysicsFreeVector.hh
        G4PhysicsLinearVector.hh
        G4PhysicsLnVector.hh
        G4PhysicsLogVector.hh
        G4PhysicsOrderedFreeVector.hh
        G4PhysicsOrderedFreeVector.icc
        G4PhysicsTable.hh
        G4PhysicsTable.icc
        G4PhysicsVector.hh
        G4PhysicsVector.icc
        G4PhysicsVectorCache.hh
        G4PhysicsVectorType.hh
        G4Physics2DVector.hh
        G4Physics2DVector.icc
        G4Physics2DVectorCache.hh
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
    SOURCES
        G4AllocatorPool.cc
        G4DataVector.cc
        G4ErrorPropagatorData.cc
        G4Exception.cc
        G4GeometryTolerance.cc
        G4LPhysicsFreeVector.cc
        G4OrderedTable.cc
        G4PhysicsFreeVector.cc
        G4PhysicsLinearVector.cc
        G4PhysicsLnVector.cc
        G4PhysicsLogVector.cc
        G4PhysicsOrderedFreeVector.cc
        G4PhysicsTable.cc
        G4PhysicsVector.cc
        G4PhysicsVectorCache.cc
        G4Physics2DVector.cc
        G4Physics2DVectorCache.cc
        G4Pow.cc
        G4ReferenceCountedHandle.cc
        G4SliceTimer.cc
        G4StateManager.cc
        G4Timer.cc
        G4UnitsTable.cc
        G4VExceptionHandler.cc
        G4VNotifier.cc
        G4VStateDependent.cc
        G4coutDestination.cc
        G4ios.cc
    GRANULAR_DEPENDENCIES
    GLOBAL_DEPENDENCIES
    LINK_LIBRARIES
        ${CLHEP_LIBRARIES}
)

# List any source specific properties here

