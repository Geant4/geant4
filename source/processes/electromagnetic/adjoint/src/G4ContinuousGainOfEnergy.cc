//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#include "G4ContinuousGainOfEnergy.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VEmModel.hh"
#include "G4VEmFluctuationModel.hh"
#include "G4VParticleChange.hh"
#include "G4UnitsTable.hh"


///////////////////////////////////////////////////////
//
G4ContinuousGainOfEnergy::G4ContinuousGainOfEnergy(const G4String& name, 
  G4ProcessType type): G4VContinuousProcess(name, type)
{

  linLossLimit=0.05;
  lossFluctuationArePossible =true;
  lossFluctuationFlag=true;
  is_integral = false;
  
}

///////////////////////////////////////////////////////
//
G4ContinuousGainOfEnergy::~G4ContinuousGainOfEnergy()
{
 
}
///////////////////////////////////////////////////////
//

void G4ContinuousGainOfEnergy::PreparePhysicsTable(
     const G4ParticleDefinition& )
{//theDirectEnergyLossProcess->PreparePhysicsTable(part);

; 
}

///////////////////////////////////////////////////////
//

void G4ContinuousGainOfEnergy::BuildPhysicsTable(const G4ParticleDefinition&)
{//theDirectEnergyLossProcess->BuildPhysicsTable(part);
;
}




///////////////////////////////////////////////////////
//
// 
G4VParticleChange* G4ContinuousGainOfEnergy::AlongStepDoIt(const G4Track& track,
                                                       const G4Step& step)
{
   
  aParticleChange.Initialize(track);
  
  // Get the actual (true) Step length
  //----------------------------------
  G4double length = step.GetStepLength();
  G4double degain  = 0.0;
 
 
 

  // Compute this for weight change after continuous energy loss
  //-------------------------------------------------------------
  G4double DEDX_before = 
  		theDirectEnergyLossProcess
				->GetDEDX(preStepKinEnergy, currentCouple);
 
 
  
  // For the fluctuation we generate a new dynamic particle with energy =preEnergy+egain
  // and then compute the fluctuation given in  the direct case.
  //-----------------------------------------------------------------------
  G4DynamicParticle* dynParticle = new G4DynamicParticle();
  *dynParticle = *(track.GetDynamicParticle());
  G4double Tkin = dynParticle->GetKineticEnergy();
  
  size_t n=1;
  if (is_integral ) n=10;
  G4double dlength= length/n; 
  for (size_t i=0;i<n;i++) {
  	G4double r = theDirectEnergyLossProcess->GetRange(Tkin, currentCouple);
   	if( dlength <= linLossLimit * r ) {
    		degain = DEDX_before*dlength;
	} 
  	else {
    		G4double x = r + length;
    		degain = theDirectEnergyLossProcess->GetKineticEnergy(x,currentCouple) - theDirectEnergyLossProcess->GetKineticEnergy(r,currentCouple);
  	}
	G4VEmModel* currentModel =  theDirectEnergyLossProcess->SelectModelForMaterial(Tkin+degain,currentMaterialIndex);
  	G4double tmax = currentModel->MaxSecondaryKinEnergy(dynParticle);
 	tmax = std::min(tmax,currentTcut);
 
  

  	// Sample fluctuations
  	//-------------------
  	
	G4double deltaE =0.;
  	if (lossFluctuationFlag ) {
      		deltaE = currentModel->GetModelOfFluctuations()->
      						SampleFluctuations(currentMaterial,dynParticle,tmax,length,degain)-degain;
  	}
	Tkin+=degain+deltaE;
	dynParticle->SetKineticEnergy(Tkin);
  	
  }
 
  

  // Corrections, which cannot be tabulated 
  // probably  this should be also changed
  // at this time it does nothing so we can leave it
  //CorrectionsAlongStep(currentCouple, dynParticle, egain, length);
  
  delete dynParticle;
 
  
  G4double DEDX_after = theDirectEnergyLossProcess->GetDEDX(Tkin, currentCouple);
  G4double weight_correction=DEDX_after/DEDX_before; //probably not needed
  weight_correction=1.;
  
  aParticleChange.ProposeEnergy(Tkin);
  
  //we still need to register in the particleChange the modification of the weight of the particle 
  G4double new_weight=weight_correction*track.GetWeight();
  aParticleChange.SetParentWeightByProcess(true);
  aParticleChange.ProposeParentWeight(new_weight);
  

  return &aParticleChange;

}
///////////////////////////////////////////////////////
//
void G4ContinuousGainOfEnergy::SetLossFluctuations(G4bool val)
{
  if(val && !lossFluctuationArePossible) return;
  lossFluctuationFlag = val;
}
