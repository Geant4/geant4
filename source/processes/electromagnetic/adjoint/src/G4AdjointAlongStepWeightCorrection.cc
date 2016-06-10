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
// $Id: G4AdjointAlongStepWeightCorrection.cc 66892 2013-01-17 10:57:59Z gunter $
//
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
 currentMaterialIndex=0;
 preStepKinEnergy=1.;
 currentCouple=0;
}

///////////////////////////////////////////////////////
//

G4AdjointAlongStepWeightCorrection::~G4AdjointAlongStepWeightCorrection()
{delete fParticleChange;
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
	
  
  

  //Caution!!!
  // It is important  to select the weight of the post_step_point
  // as the current weight and not the weight of the track, as t
  // the  weight of the track is changed after having applied all
  // the along_step_do_it.

  // G4double new_weight=weight_correction*track.GetWeight(); //old
  G4double new_weight=weight_correction*step.GetPostStepPoint()->GetWeight();



  //if (weight_correction >2.) new_weight=1.e-300;
  
  
  //The following test check for zero weight.
  //This happens after weight correction of gamma for photo electric effect.
  //When the new weight is 0 it will be later on consider as nan by G4.
  //Therefore we do put a lower limit of 1.e-300. for new_weight 
  //Correction by L.Desorgher on 15 July 2009 
  if (new_weight==0 || (new_weight<=0 && new_weight>0)){
		//G4cout<<new_weight<<'\t'<<weight_correction<<'\t'<<track.GetWeight()<<G4endl;
		new_weight=1.e-300;
  }
  
  //G4cout<<new_weight<<'\t'<<weight_correction<<'\t'<<track.GetWeight()<<G4endl;
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
