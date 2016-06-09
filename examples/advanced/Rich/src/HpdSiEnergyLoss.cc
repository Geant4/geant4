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
// Rich advanced example for Geant4
// HpdSiEnergyLoss.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include "HpdSiEnergyLoss.hh"
#include "G4Material.hh"
#include "Randomize.hh"
#include "RichTbMaterialParameters.hh"
#include "RichTbAnalysisManager.hh"

HpdSiEnergyLoss::HpdSiEnergyLoss(const G4String& materialName,
                                 const G4String& processName)
  : G4VEnergyLoss(processName),
    MinKineticEnergy(1.*eV),MipEnergy(30000.0*eV),
    finalRangeforStep(0.15*mm)  { 
  ElossMaterialName= materialName;
  
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  G4int numberOfMat = G4Material::GetNumberOfMaterials();
  
  G4int iMat=0;
  for(iMat=0;iMat<numberOfMat;iMat++) {
    if ( materialName == (*theMaterialTable)[iMat]->GetName()){
      fMatIndex=(*theMaterialTable)[iMat]->GetIndex();
      break;
    }
    
  }
  if(iMat >= numberOfMat ) {
    G4Exception("Invalid material Name in HpdSiEnergyLoss constructor" );
  }
}

HpdSiEnergyLoss::~HpdSiEnergyLoss() { }

G4bool HpdSiEnergyLoss::IsApplicable(const G4ParticleDefinition& 
				     aParticleType) {
   return(aParticleType.GetPDGCharge()!= 0.);
}

G4double HpdSiEnergyLoss::GetContinuousStepLimit(const G4Track& track,
                                 G4double,
                                 G4double,
                                 G4double& ){

  G4double  RangeForStep =  finalRangeforStep;

  if( fMatIndex != G4int(track.GetMaterial() -> GetIndex()) ) { 
    RangeForStep = DBL_MAX;
  }

  return RangeForStep;
}

G4double  HpdSiEnergyLoss::GetMeanFreePath(const G4Track&,
                         G4double,
                         G4ForceCondition* condition) {
  // return infinity so that it does nothing.
  *condition = NotForced;
  return DBL_MAX;

}

G4VParticleChange*  HpdSiEnergyLoss::PostStepDoIt(const G4Track& aTrack,
                                 const G4Step& aStep) {
  // Do nothing
   aParticleChange.Initialize(aTrack) ;
  return G4VContinuousDiscreteProcess::PostStepDoIt(aTrack,aStep);
   
}

G4VParticleChange* HpdSiEnergyLoss::AlongStepDoIt(const G4Track& aTrack,
				 const G4Step&) {


#ifdef G4ANALYSIS_USE
  RichTbAnalysisManager * analysis = RichTbAnalysisManager::getInstance();
#endif


  aParticleChange.Initialize(aTrack);
  G4int aMaterialIndex = aTrack.GetMaterial()->GetIndex();
  if(fMatIndex != aMaterialIndex ) {
    return &aParticleChange;
  }

  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  G4double aKinEnergyInit = aParticle->GetKineticEnergy();
  G4double Eloss, aKinEnergyFinal;
  if(aKinEnergyInit < MinKineticEnergy ) {  Eloss=0.0 ; }
  else if( aKinEnergyInit < MipEnergy ) {Eloss= aKinEnergyInit ;}
  else { Eloss = MipEnergy; }
 
  aKinEnergyFinal=aKinEnergyInit-Eloss;


  //In the G4example the backscattering is implemented in
  //an adhoc manner as done  below. It simply causes an
  // efficiency loss.

  G4double bckratio;
  if( SignalToNoiseInData > 0.0 ){
    bckratio =  NsigmaInPedCut/ SignalToNoiseInData ;

  }
  G4double Effs = 1.0 -  backscaprob * bckratio;
  G4double Randbsk =  G4UniformRand();
  if(Randbsk <= Effs && Eloss > 0.0 ) {
           aParticleChange.ProposeLocalEnergyDeposit(Eloss);

#ifdef G4ANALYSIS_USE
  
	   analysis->bumpNumHitInSi();

#endif

    // G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
    // G4String tpreVol = pPreStepPoint -> GetPhysicalVolume()->GetName();
    // G4int tpreVP =  pPreStepPoint -> GetPhysicalVolume()->GetCopyNo();
  }
  if (aKinEnergyFinal <= MinKineticEnergy ) {
       aParticleChange.ProposeTrackStatus(fStopAndKill);
  
   }else { 
     aParticleChange.ProposeEnergy(aKinEnergyFinal);

   }
  return &aParticleChange;

}








