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
// Authors: Susanna Guatelli, susanna@uow.edu.au,
// Authors: Jeremy Davis, jad028@uowmail.edu.au
//

#include "G4ios.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "SteppingAction.hh"
#include "AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

SteppingAction::SteppingAction(AnalysisManager* pAnalysis)
{ 
analysis = pAnalysis;
fSecondary = 0; 
}

SteppingAction::~SteppingAction()
{ 
}

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
  G4SteppingManager*  steppingManager = fpSteppingManager;
  G4Track* theTrack = aStep -> GetTrack();

  // check if it is alive
  if(theTrack-> GetTrackStatus() == fAlive) { return; }

  // Retrieve the secondary particles
    fSecondary = steppingManager -> GetfSecondary();

#ifdef ANALYSIS_USE  
   for(size_t lp1=0;lp1<(*fSecondary).size(); lp1++)
     { 
       // Retrieve the info about the generation of secondary particles 
       G4String volumeName = (*fSecondary)[lp1] -> GetVolume() -> GetName(); // volume where secondary was generated 
       G4String secondaryParticleName =  (*fSecondary)[lp1]->GetDefinition() -> GetParticleName();  // name of the secondary
       G4double secondaryParticleKineticEnergy =  (*fSecondary)[lp1] -> GetKineticEnergy(); // kinetic energy
      // G4String process = (*fSecondary)[lp1]-> GetCreatorProcess()-> GetProcessName();   // process creating it
       G4double charge = (*fSecondary)[lp1] -> GetDynamicParticle() -> GetDefinition() -> GetPDGCharge();
       G4int AA = (*fSecondary)[lp1] -> GetDynamicParticle() -> GetDefinition() -> GetBaryonNumber();
          
     if (volumeName == "SV_phys1") 
//Testing with larger volume!
//if (volumeName == "DiaVol_phys")
	 {
	   if ((secondaryParticleName == "proton") ||
               (secondaryParticleName == "neutron")||
               (secondaryParticleName == "alpha") ||
               (secondaryParticleName == "deuton") || 
               (secondaryParticleName == "triton") || 
               (secondaryParticleName == "He3") || 
	       (secondaryParticleName =="GenericIon"))
               analysis -> FillSecondaries(AA, charge, secondaryParticleKineticEnergy/MeV); 	
          }
   }
#endif	 
}

