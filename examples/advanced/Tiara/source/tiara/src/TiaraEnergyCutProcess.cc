
#include "TiaraEnergyCutProcess.hh"

TiaraEnergyCutProcess::TiaraEnergyCutProcess(G4double minEnergyCut)
 : 
  G4VProcess("TiaraEnergyCutProcess"), 
  fMinEnergyCut(minEnergyCut)
{
  G4VProcess::pParticleChange = new G4ParticleChange;
  if (!G4VProcess::pParticleChange) {
    G4Exception("ERROR:TiaraEnergyCutProcess::TiaraEnergyCutProcess new failed to create G4ParticleChange!");
  }
}

TiaraEnergyCutProcess::~TiaraEnergyCutProcess()
{
  delete pParticleChange;
}

G4double TiaraEnergyCutProcess::
PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				     G4double   previousStepSize,
				     G4ForceCondition* condition)
{
  *condition = Forced;
  return kInfinity;
}
  
G4VParticleChange * 
TiaraEnergyCutProcess::PostStepDoIt(const G4Track& aTrack, const G4Step &aStep)
{
  pParticleChange->Initialize(aTrack);

  if (aTrack.GetKineticEnergy() < fMinEnergyCut) {
    pParticleChange->SetStatusChange(fStopAndKill);
  }
  
  return G4VProcess::pParticleChange;
}

const G4String &TiaraEnergyCutProcess::GetName() const {
  return theProcessName;
}


G4double TiaraEnergyCutProcess::
AlongStepGetPhysicalInteractionLength(const G4Track&,
				      G4double  ,
				      G4double  ,
				      G4double& ,
				      G4GPILSelection*) {
  return -1.0;
}

G4double TiaraEnergyCutProcess::
AtRestGetPhysicalInteractionLength(const G4Track&,
				   G4ForceCondition*) {
  return -1.0;
}

G4VParticleChange* TiaraEnergyCutProcess::AtRestDoIt(const G4Track&,
					       const G4Step&) {
  return 0;
}

G4VParticleChange* TiaraEnergyCutProcess::AlongStepDoIt(const G4Track&,
						  const G4Step&) {
  return 0;
}

