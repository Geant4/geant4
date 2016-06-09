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
