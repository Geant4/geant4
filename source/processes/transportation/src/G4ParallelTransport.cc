#include "G4ParallelTransport.hh"
#include "G4VPGeoDriver.hh"
#include "G4VParallelStepper.hh"
#include "G4Pstring.hh"




G4ParallelTransport::G4ParallelTransport(G4VPGeoDriver &pgeodriver,
					 G4VParallelStepper &aStepper,
					 const G4String &aName) : 
  G4VProcess(aName), 
  fPgeodriver(pgeodriver),
  fPStepper(aStepper),
  fCrossBoundary(false){
  fParticleChange = new G4ParticleChange;
  G4VProcess::pParticleChange = fParticleChange;
}
G4ParallelTransport::~G4ParallelTransport(){
  delete fParticleChange;
}
  
G4double 
G4ParallelTransport::
PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				     G4double   previousStepSize,
				     G4ForceCondition* condition){

  if (aTrack.GetCurrentStepNumber()!=1 &&  previousStepSize==0.) {
    G4String m = "GPIL CurrentStepNumber() = ";
    m += str(aTrack.GetCurrentStepNumber()); 
    m += " &&  previousStepSize==0\n";
    m +=  "pos: " + str(aTrack.GetPosition()); 
    m +=  ", dir: " + str(aTrack.GetMomentumDirection());
    Warning(m);
  }

  G4double stepLength;

  // if this function is called on a new track let the navigator know about it
  G4bool initStep = (aTrack.GetCurrentStepNumber() == 1);

  if (initStep) {
    stepLength = fPgeodriver.
      ComputeStepLengthInit(aTrack.GetPosition(),
			    aTrack.GetMomentumDirection());
    fPStepper.Init(fPgeodriver.GetCurrentTouchableKey());
    fCrossBoundary = false;
  }
  else if (fCrossBoundary) {
    stepLength = fPgeodriver.
      ComputeStepLengthCrossBoundary(aTrack.GetPosition(),
				     aTrack.GetMomentumDirection());
    fPStepper.UnSetCrossBoundary();
    fCrossBoundary = false;
  }
  else {
    stepLength = fPgeodriver.
      ComputeStepLengthInVolume(aTrack.GetPosition(),
				aTrack.GetMomentumDirection());
  }

    
  
  *condition = NotForced;
  return stepLength;
}

G4VParticleChange * 
G4ParallelTransport::PostStepDoIt(const G4Track& aTrack,
				  const G4Step& aStep){
  if (aStep.GetStepLength() == 0.) {
    G4String m = "G4PArallelTransport::PostStepDoIt: StepLength() == 0.\n";
    m += "pos: " + str(aTrack.GetPosition()) + ", " 
      + "dir: " +  str(aTrack.GetMomentumDirection()) + "\n";
    Warning(m);
  }

  fParticleChange->Initialize(aTrack);

  fCrossBoundary = true;
  fPgeodriver.LocateOnBoundary(aTrack.GetPosition(), 
			       aTrack.GetMomentumDirection());
  fPStepper.Update(fPgeodriver.GetCurrentTouchableKey());

  return fParticleChange;
}
    

void G4ParallelTransport::Warning(const G4String &m){
  G4cout << "WARNING: G4ParallelTransport:: " << G4endl;
  G4cout << m << G4endl;
}
