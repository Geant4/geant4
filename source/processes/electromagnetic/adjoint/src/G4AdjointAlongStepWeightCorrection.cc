#include "G4AdjointAlongStepWeightCorrection.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VParticleChange.hh"
#include "G4AdjointCSManager.hh"


///////////////////////////////////////////////////////
//

G4AdjointAlongStepWeightCorrection::G4AdjointAlongStepWeightCorrection(const G4String& name, 
  G4ProcessType type): G4VContinuousProcess(name, type)
{fParticleChange = new G4ParticleChange();
}

///////////////////////////////////////////////////////
//

G4AdjointAlongStepWeightCorrection::~G4AdjointAlongStepWeightCorrection()
{; 
}


///////////////////////////////////////////////////////
//

void G4AdjointAlongStepWeightCorrection::PreparePhysicsTable(
     const G4ParticleDefinition& )
{
; 

}

///////////////////////////////////////////////////////
//

void G4AdjointAlongStepWeightCorrection::BuildPhysicsTable(const G4ParticleDefinition& )
{;
}




///////////////////////////////////////////////////////
//
G4VParticleChange* G4AdjointAlongStepWeightCorrection::AlongStepDoIt(const G4Track& track,
                                                       const G4Step& step)
{
   
  fParticleChange->Initialize(track);
  
  // Get the actual (true) Step length
  //----------------------------------
  G4double length = step.GetStepLength();
 

  G4double Tkin = step.GetPostStepPoint()->GetKineticEnergy();
  G4ParticleDefinition* thePartDef= const_cast<G4ParticleDefinition*>  (track.GetDynamicParticle()->GetDefinition());
  G4double weight_correction=G4AdjointCSManager::GetAdjointCSManager()->GetContinuousWeightCorrection(thePartDef,
  									preStepKinEnergy,Tkin, currentCouple,length);
	
 
  
  G4double new_weight=weight_correction*track.GetWeight();
  fParticleChange->SetParentWeightByProcess(false);
  fParticleChange->SetSecondaryWeightByProcess(false);
  fParticleChange->ProposeParentWeight(new_weight);
  

  return fParticleChange;

}
