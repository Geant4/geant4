#include "G4MassImportanceProcess.hh"
#include "G4VImportanceAlgorithm.hh"
#include "G4ImportanceFinder.hh"
#include "G4PTouchableKey.hh"

G4MassImportanceProcess::
G4MassImportanceProcess(const G4VImportanceAlgorithm &aImportanceAlgorithm,
			const G4VIStore &aIstore,
			const G4String &aName):
  fImportanceAlgorithm(aImportanceAlgorithm),
  fImportanceFinder(new G4ImportanceFinder(aIstore)){
  fParticleChange = new G4ParticleChange;
  G4VProcess::pParticleChange = fParticleChange;
}

G4MassImportanceProcess::~G4MassImportanceProcess() {
  delete fImportanceFinder;
  delete fParticleChange;
}

G4double 
G4MassImportanceProcess::
PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				     G4double   previousStepSize,
				     G4ForceCondition* condition){
  *condition = Forced;
  return kInfinity;
}
  
G4VParticleChange * 
G4MassImportanceProcess::PostStepDoIt(const G4Track &aTrack,
				      const G4Step &aStep){
  fParticleChange->Initialize(aTrack);
  if (aStep.GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
    if (aTrack.GetTrackStatus()==fStopAndKill) {
      G4cout << "G4MassImportanceProcess::PostStepDoIt StopAndKill" << G4endl;
    }

    G4StepPoint *prepoint = aStep.GetPreStepPoint();
    G4StepPoint *postpoint = aStep.GetPostStepPoint();
  

    G4PTouchableKey prekey(*(prepoint->GetPhysicalVolume()), 
			 prepoint->GetTouchable()->GetReplicaNumber());
    G4PTouchableKey postkey(*(postpoint->GetPhysicalVolume()), 
			  postpoint->GetTouchable()->GetReplicaNumber());

    G4Nsplit_Weight nw = fImportanceAlgorithm.
      Calculate(fImportanceFinder->
		GetIPre_over_IPost(prekey, postkey),
		aTrack.GetWeight());
    fImportancePostStepDoIt.DoIt(aTrack, fParticleChange, nw);
  }
  return fParticleChange;
}
