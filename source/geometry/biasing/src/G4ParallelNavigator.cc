#include "G4ParallelNavigator.hh"

#include "G4VTouchable.hh"
#include "G4Navigator.hh"
#include "G4PTouchableKey.hh"
#include "G4VParallelStepper.hh"
#include "g4std/strstream"

using namespace std;

G4ParallelNavigator::G4ParallelNavigator(G4VPhysicalVolume &aWorldVolume): 
  fNavigator(*(new G4Navigator)),
  fNlocated(0){
  fNavigator.SetWorldVolume(&aWorldVolume);
  fCurrentTouchable = fNavigator.CreateTouchableHistory();
}
G4ParallelNavigator::~G4ParallelNavigator(){
  delete fCurrentTouchable;
  delete &fNavigator;
}

// public functions

G4PTouchableKey G4ParallelNavigator::
LocateOnBoundary(const G4ThreeVector &aPosition, 
		 const G4ThreeVector &aDirection) {
  Locate(aPosition, aDirection, true);
  // since the track crosses a boundary ipdate stepper 
  return GetCurrentTouchableKey();
}

G4double G4ParallelNavigator::
ComputeStepLengthInit(const G4ThreeVector &aPosition, 
		      const G4ThreeVector &aDirection){
  G4double newSafety;
  G4double stepLength;

  // initialization
  fNlocated = 0;
  
  // first try
  Locate(aPosition, aDirection, false);
  stepLength  = fNavigator.ComputeStep( aPosition, aDirection,
					kInfinity, newSafety);
  if (stepLength<=0.0) {
    // second try with schifted position
    stepLength  = ComputeStepLengthShifted("ComputeStepLengthInit",
					   aPosition, aDirection);
  } 
  return stepLength;
}

G4double G4ParallelNavigator::
ComputeStepLengthCrossBoundary(const G4ThreeVector &aPosition, 
			       const G4ThreeVector &aDirection){
  if (fNlocated == 0 ) {
    Error("ComputeStepLengthCrossBoundary: no location done before call",
	  aPosition, aDirection);
  }

  // if the track is on boundary the location was done
  // in the DOIT of the ImportanceSampler

  G4double newSafety;    
  G4double stepLength = fNavigator.ComputeStep( aPosition, aDirection,
				       kInfinity, newSafety);
  if (stepLength<=0.) {
    // second try with schifted position
    G4cout << "Warning: G4ParallelNavigator::ComputeStepLengthCrossBoundary: "
	   << G4endl;
    G4cout << "stepLength<=0." << ", pos = " << aPosition 
	   << ", dir = " << aDirection << G4endl;
    G4cout << "try ComputeStepLengthShifted" << G4endl;
    stepLength  = ComputeStepLengthShifted("ComputeStepLengthCrossBoundary",
					   aPosition, aDirection);
  }
  return stepLength;
}

G4double G4ParallelNavigator::
ComputeStepLengthInVolume(const G4ThreeVector &aPosition, 
			  const G4ThreeVector &aDirection){
  if (fNlocated == 0 ) {
    Error("ComputeStepLengthInVolumeno location done before call", 
	  aPosition, aDirection);
  }

  // if the track is not on the boundary and it's
  // not the first step, the location must be inside the 
  // volume 
  fNavigator.LocateGlobalPointWithinVolume(aPosition);
  G4double newSafety;
  G4double stepLength = fNavigator.ComputeStep( aPosition, aDirection,
				       kInfinity, newSafety);
  if (stepLength<=0.) {
    // second try with schifted position
    G4cout << "Warning: G4ParallelNavigator::ComputeStepLengthInVolume: "
	   << G4endl;
    G4cout << "stepLength<=0." << ", pos = " << aPosition 
	   << ", dir = " << aDirection << G4endl;
    G4cout << "try ComputeStepLengthShifted" << G4endl;
    stepLength  = ComputeStepLengthShifted("ComputeStepLengthInVolume",
					   aPosition, aDirection);
  }
  return stepLength;
}

// private functions

void G4ParallelNavigator::Locate(const G4ThreeVector &aPosition, 
				   const G4ThreeVector &aDirection,
				   G4bool histsearch){
  fNavigator.SetGeometricallyLimitedStep();
  fNavigator.
    LocateGlobalPointAndUpdateTouchable( aPosition, aDirection,
                                         fCurrentTouchable, histsearch);
  fNlocated++;
  
  return;
  
}



G4double G4ParallelNavigator::
ComputeStepLengthShifted(const G4String &m,
			 const G4ThreeVector &aPosition, 
			 const G4ThreeVector &aDirection){
  G4ThreeVector shift_pos(aPosition);
  shift_pos+=G4ThreeVector(Shift(aDirection.x()), 
			   Shift(aDirection.y()), 
			   Shift(aDirection.z()));
  Locate(shift_pos, aDirection, false);
  G4double newSafety;
  G4double stepLength  = fNavigator.ComputeStep( shift_pos, aDirection,
					kInfinity, newSafety);
  if (stepLength<=0.0) {
    ostrstream os;
    os << "ComputeStepLengthShifted: called by " << m << "\n";
    os << "shifted positio in second try: " << shift_pos << '\0' << G4endl;
    Error(os.str() , aPosition, aDirection);
  } 
  return stepLength;
}


G4PTouchableKey G4ParallelNavigator::GetCurrentTouchableKey() const {
  return G4PTouchableKey(*fCurrentTouchable->GetVolume(),
			 fCurrentTouchable->GetReplicaNumber());
}

void G4ParallelNavigator::
Error(const G4String &m, const G4ThreeVector &pos, const G4ThreeVector &dir){
  G4cout << "ERROE: in G4ParallelNavigator::" <<  m << G4endl;
  G4cout << "aPosition: " << pos << G4endl;
  G4cout << "dir: " << dir << G4endl;
  exit(1);
}





