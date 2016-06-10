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
// $Id: G4ContinuousGainOfEnergy.cc 91870 2015-08-07 15:21:40Z gcosmo $
//

#include "G4ContinuousGainOfEnergy.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VEmModel.hh"
#include "G4VEmFluctuationModel.hh"
#include "G4VParticleChange.hh"
#include "G4AdjointCSManager.hh"
#include "G4LossTableManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

///////////////////////////////////////////////////////
//
G4ContinuousGainOfEnergy::G4ContinuousGainOfEnergy(const G4String& name, 
  G4ProcessType type): G4VContinuousProcess(name, type)
{


  linLossLimit=0.05;
  lossFluctuationArePossible =true;
  lossFluctuationFlag=true;
  is_integral = false;
  
  //Will be properly set in SetDirectParticle()
  IsIon=false;
  massRatio =1.;
  chargeSqRatio=1.;
  preStepChargeSqRatio=1.;
  
  //Some initialization
  currentCoupleIndex=9999999;
  currentCutInRange=0.;
  currentMaterialIndex=9999999;
  currentTcut=0.;
  preStepKinEnergy=0.;
  preStepRange=0.;
  preStepScaledKinEnergy=0.;
  
  currentCouple=0;  
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
void  G4ContinuousGainOfEnergy::SetDirectParticle(G4ParticleDefinition* p)
{theDirectPartDef=p;
 if (theDirectPartDef->GetParticleType()== "nucleus") {
 	 IsIon=true;
	 massRatio = proton_mass_c2/theDirectPartDef->GetPDGMass();
	 G4double q=theDirectPartDef->GetPDGCharge();
	 chargeSqRatio=q*q;
	
	 
 }
 
}

///////////////////////////////////////////////////////
//
// 
G4VParticleChange* G4ContinuousGainOfEnergy::AlongStepDoIt(const G4Track& track,
                                                       const G4Step& step)
{
   
  //Caution in this method the  step length should be the true step length
  // A problem is that this is compute by the multiple scattering that does not know the energy at the end of the adjoint step. This energy is used during the 
  //Forward sim. Nothing we can really do against that at this time. This is inherent to the MS method
  //
  
  
  
  aParticleChange.Initialize(track);
  
  // Get the actual (true) Step length
  //----------------------------------
  G4double length = step.GetStepLength();
  G4double degain  = 0.0;
  
  
 
  // Compute this for weight change after continuous energy loss
  //-------------------------------------------------------------
  G4double DEDX_before = theDirectEnergyLossProcess->GetDEDX(preStepKinEnergy, currentCouple);
   
  
  
  // For the fluctuation we generate a new dynamic particle with energy =preEnergy+egain
  // and then compute the fluctuation given in  the direct case.
  //-----------------------------------------------------------------------
  G4DynamicParticle* dynParticle = new G4DynamicParticle();
  *dynParticle = *(track.GetDynamicParticle());
  dynParticle->SetDefinition(theDirectPartDef);
  G4double Tkin = dynParticle->GetKineticEnergy(); 


  size_t n=1;
  if (is_integral ) n=10;
  n=1;
  G4double dlength= length/n; 
  for (size_t i=0;i<n;i++) {
  	if (Tkin != preStepKinEnergy && IsIon) {
  		chargeSqRatio =  currentModel->GetChargeSquareRatio(theDirectPartDef,currentMaterial,Tkin);
		theDirectEnergyLossProcess->SetDynamicMassCharge(massRatio,chargeSqRatio); 
	
  	}
  
  	G4double r = theDirectEnergyLossProcess->GetRange(Tkin, currentCouple);
   	if( dlength <= linLossLimit * r ) {
    		degain = DEDX_before*dlength;		
	} 
  	else {
    		G4double x = r + dlength;
    		//degain = theDirectEnergyLossProcess->GetKineticEnergy(x,currentCouple) - theDirectEnergyLossProcess->GetKineticEnergy(r,currentCouple);
		G4double E = theDirectEnergyLossProcess->GetKineticEnergy(x,currentCouple);
		if (IsIon){
			chargeSqRatio =  currentModel->GetChargeSquareRatio(theDirectPartDef,currentMaterial,E);
			theDirectEnergyLossProcess->SetDynamicMassCharge(massRatio,chargeSqRatio);
			G4double x1= theDirectEnergyLossProcess->GetRange(E, currentCouple);

                        // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
                        G4int ii=0;
                        const G4int iimax = 100;
			while (std::abs(x-x1)>0.01*x) {
				E = theDirectEnergyLossProcess->GetKineticEnergy(x,currentCouple);
				chargeSqRatio =  currentModel->GetChargeSquareRatio(theDirectPartDef,currentMaterial,E);
				theDirectEnergyLossProcess->SetDynamicMassCharge(massRatio,chargeSqRatio);
				x1= theDirectEnergyLossProcess->GetRange(E, currentCouple);
				++ii;
			        if(ii >= iimax) { break; }
			} 
		}
		
		degain=E-Tkin;	
		
		
		
  	}
	//G4cout<<degain<<G4endl;
  	G4double tmax = currentModel->MaxSecondaryKinEnergy(dynParticle);
 	tmax = std::min(tmax,currentTcut);
 	
	
	dynParticle->SetKineticEnergy(Tkin+degain);

	// Corrections, which cannot be tabulated for ions
	//----------------------------------------
	G4double esecdep=0;//not used in most models
	currentModel->CorrectionsAlongStep(currentCouple, dynParticle, degain,esecdep, dlength); 

  	// Sample fluctuations
  	//-------------------
	
  	
	G4double deltaE =0.;
  	if (lossFluctuationFlag ) {
      	  deltaE = currentModel->GetModelOfFluctuations()->
      	    SampleFluctuations(currentCouple,dynParticle,tmax,dlength,degain)-degain;
  	}
	
	G4double egain=degain+deltaE;
	if (egain <=0) egain=degain;
	Tkin+=egain;
	dynParticle->SetKineticEnergy(Tkin);
 }
 
 
  

  
  delete dynParticle;
 
  if (IsIon){
	chargeSqRatio =  currentModel->GetChargeSquareRatio(theDirectPartDef,currentMaterial,Tkin);
	theDirectEnergyLossProcess->SetDynamicMassCharge(massRatio,chargeSqRatio);
		
  }
  
  G4double DEDX_after = theDirectEnergyLossProcess->GetDEDX(Tkin, currentCouple);
  
  
  G4double weight_correction=DEDX_after/DEDX_before;
 
  
  aParticleChange.ProposeEnergy(Tkin);


  //Caution!!!
  // It is important  to select the weight of the post_step_point
  // as the current weight and not the weight of the track, as t
  // the  weight of the track is changed after having applied all
  // the along_step_do_it.

  // G4double new_weight=weight_correction*track.GetWeight(); //old
  G4double new_weight=weight_correction*step.GetPostStepPoint()->GetWeight();
  aParticleChange.SetParentWeightByProcess(false);
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
///////////////////////////////////////////////////////
//



G4double G4ContinuousGainOfEnergy::GetContinuousStepLimit(const G4Track& track,
                G4double , G4double , G4double& )
{ 
  G4double x = DBL_MAX;
  x=.1*mm;
 
 
  DefineMaterial(track.GetMaterialCutsCouple());
 
  preStepKinEnergy = track.GetKineticEnergy(); 
  preStepScaledKinEnergy = track.GetKineticEnergy()*massRatio;
  currentModel = theDirectEnergyLossProcess->SelectModelForMaterial(preStepScaledKinEnergy,currentCoupleIndex);
  G4double emax_model=currentModel->HighEnergyLimit();
  if (IsIon) {
  	chargeSqRatio =  currentModel->GetChargeSquareRatio(theDirectPartDef,currentMaterial,preStepKinEnergy);
	preStepChargeSqRatio = chargeSqRatio;
	theDirectEnergyLossProcess->SetDynamicMassCharge(massRatio,preStepChargeSqRatio);
  } 
  
  
  G4double maxE =1.1*preStepKinEnergy;
  /*if (preStepKinEnergy< 0.05*MeV) maxE =2.*preStepKinEnergy;
  else if (preStepKinEnergy< 0.1*MeV) maxE =1.5*preStepKinEnergy;
  else if (preStepKinEnergy< 0.5*MeV) maxE =1.25*preStepKinEnergy;*/
   
  if (preStepKinEnergy < currentTcut) maxE = std::min(currentTcut,maxE);
 
  maxE=std::min(emax_model*1.001,maxE);
  	
  preStepRange = theDirectEnergyLossProcess->GetRange(preStepKinEnergy, currentCouple);
  
  if (IsIon) {
  	G4double chargeSqRatioAtEmax = currentModel->GetChargeSquareRatio(theDirectPartDef,currentMaterial,maxE);
	theDirectEnergyLossProcess->SetDynamicMassCharge(massRatio,chargeSqRatioAtEmax);
  }	
  
  G4double r1 = theDirectEnergyLossProcess->GetRange(maxE, currentCouple);
  
  if (IsIon) theDirectEnergyLossProcess->SetDynamicMassCharge(massRatio,preStepChargeSqRatio);
  
  

  x=r1-preStepRange;
  x=std::max(r1-preStepRange,0.001*mm);
 
  return x;
  
 
}
#include "G4EmCorrections.hh"
///////////////////////////////////////////////////////
//

void G4ContinuousGainOfEnergy::SetDynamicMassCharge(const G4Track& ,G4double energy)
{ 

  G4double ChargeSqRatio= G4LossTableManager::Instance()->EmCorrections()->EffectiveChargeSquareRatio(theDirectPartDef,currentMaterial,energy); 
  if (theDirectEnergyLossProcess) theDirectEnergyLossProcess->SetDynamicMassCharge(massRatio,ChargeSqRatio);
}
