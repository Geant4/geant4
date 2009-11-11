#include "G4AdjointAlongStepWeightCorrection.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VParticleChange.hh"
#include "G4AdjointCSManager.hh"

#ifdef WIN32
#include <float.h>
#endif

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
  
  //if (weight_correction >2.) new_weight=1.e-300;
  
  
  //The following test check for zero weight.
  //This happens after weight correction of gamma for photo electric effect.
  //When the new weight is 0 it will be later on consider as nan by G4.
  //Therefore we do put a lower limit of 1.e-300. for new_weight 
  //Correction by L.Desorgher on 15 July 2009 
#ifdef WIN32
  if (!!_isnan(new_weight) || new_weight==0){
#else
  if (std::isnan(new_weight) || new_weight==0){
#endif
		//std::cout<<new_weight<<'\t'<<weight_correction<<'\t'<<track.GetWeight()<<std::endl;
		new_weight=1.e-300;
  }
  
  //std::cout<<new_weight<<'\t'<<weight_correction<<'\t'<<track.GetWeight()<<std::endl;
  fParticleChange->SetParentWeightByProcess(false);
  fParticleChange->SetSecondaryWeightByProcess(false);
  fParticleChange->ProposeParentWeight(new_weight);
  

  return fParticleChange;

}
///////////////////////////////////////////////////////
//
G4double G4AdjointAlongStepWeightCorrection::GetContinuousStepLimit(const G4Track& track,
                G4double , G4double , G4double& )
{ 
  G4double x = DBL_MAX;
  DefineMaterial(track.GetMaterialCutsCouple());
  preStepKinEnergy = track.GetKineticEnergy();
  return x;
}
