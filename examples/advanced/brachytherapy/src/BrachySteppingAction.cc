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
//
//
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//

#include "G4ios.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"

#ifdef ANALYSIS_USE
#include "BrachyAnalysisManager.hh"
#endif

#include "BrachySteppingAction.hh"
#include "G4SystemOfUnits.hh"

BrachySteppingAction::BrachySteppingAction()
{ 

}

BrachySteppingAction::~BrachySteppingAction()
{ 
}

void BrachySteppingAction::UserSteppingAction(const G4Step* aStep)
{

// Retrieve the spectrum of gamma emitted in the Radioactive Decay
// and store it in a 1D histogram

  G4SteppingManager*  steppingManager = fpSteppingManager;
  G4Track* theTrack = aStep-> GetTrack();

  // check if it is alive
  if(theTrack-> GetTrackStatus() == fAlive) {return;}
   
  // G4cout << "Start secondariessss" << G4endl;  
  // Retrieve the secondary particles
  G4TrackVector* fSecondary = steppingManager -> GetfSecondary();
     
   for(size_t lp1=0;lp1<(*fSecondary).size(); lp1++)
   { 
     // Retrieve particle
     const G4ParticleDefinition* particleName = (*fSecondary)[lp1] -> GetDefinition();     

     if (particleName == G4Gamma::Definition())
     {
      G4String process = (*fSecondary)[lp1]-> GetCreatorProcess()-> GetProcessName();  
      
      // Retrieve the process originating it
      // G4cout << "creator process " << process << G4endl;
        if (process == "RadioactiveDecay")
         {
#ifdef ANALYSIS_USE  
          BrachyAnalysisManager* analysis = BrachyAnalysisManager::GetInstance();
          G4double energy = (*fSecondary)[lp1]  -> GetKineticEnergy();
          // Store the initial energy of particles in a 1D histogram
          analysis -> FillPrimaryParticleHistogram(energy/keV);
#endif
         }
      }
   } 
}
