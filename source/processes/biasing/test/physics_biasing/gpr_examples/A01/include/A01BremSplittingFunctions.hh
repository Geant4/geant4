#ifndef A01BREMSPLITTINGFUNCTIONS_HH
#define A01BREMSPLITTINGFUNCTIONS_HH
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
// $Id: A01BremSplittingFunctions.hh,v 1.1 2007-09-17 08:48:06 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Jane Tinslay, May 2007
//
#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include <assert.h>
#include <vector>

namespace A01BremSplittingFunctions {

  G4VParticleChange* BremSplitting(G4GPRProcessWrappers::G4GPRDiscreteDoIt& original, 
				   const G4Track& track, const G4Step& step)
  {
    G4cout<<"jane doing brem splitting for "<<original.GetIdentifier()<<" "<<track.GetTrackID()<<G4endl;

    unsigned nSplit = 100;
    
    G4VParticleChange* particleChange = original(track, step);
    assert (0 != particleChange);
    
    unsigned i(0);
    G4double weight = track.GetWeight()/nSplit;
    
  // Secondary store
    std::vector<G4Track*> secondaries;
    secondaries.reserve(nSplit);
    
    // Loop over PostStepDoIt method to generate multiple secondaries.
    for (i=0; i<nSplit; i++) {    
      particleChange = original(track, step);
      
      assert (0 != particleChange);
      particleChange->SetVerboseLevel(0);
      
      G4int j(0);
      
      for (j=0; j<particleChange->GetNumberOfSecondaries(); j++) {
	secondaries.push_back(new G4Track(*(particleChange->GetSecondary(j))));
      }
    }	
    
    // Configure particleChange to handle multiple secondaries. Other 
    // data is unchanged
    particleChange->SetNumberOfSecondaries(secondaries.size());
    particleChange->SetSecondaryWeightByProcess(true);
    
    // Add all secondaries 
    std::vector<G4Track*>::iterator iter = secondaries.begin();
    
    while (iter != secondaries.end()) {
      G4Track* myTrack = *iter;
      myTrack->SetWeight(weight);
      
      // particleChange takes ownership
      particleChange->AddSecondary(myTrack); 
      
      iter++;
    }
    
    return particleChange;
  }

  G4VParticleChange* Roulette(G4GPRProcessWrappers::G4GPRDiscreteDoIt& original, 
			      const G4Track& track, const G4Step& step)
  {
    G4cout<<"jane doing roulette for "<<original.GetIdentifier()<<" for track "<<track.GetTrackID()<<G4endl;

    // Just do regular processing if brem splitting is not activated
    G4VParticleChange* particleChange = original(track, step);
    assert (0 != particleChange);
    
    G4int nSecondaries = particleChange->GetNumberOfSecondaries();
    
    if (0 == nSecondaries) return particleChange;
    
    particleChange->SetVerboseLevel(0);
    
    // Vector of surviving secondaries
    std::vector<G4Track*> survivors;
    
    // Russian roulette probability threshold
    G4int fFactor = 100;
    G4double inverseFactor = 1.0/static_cast<G4double>(fFactor);

    G4int i(0);
    
    for (i = 0; i<nSecondaries; i++) {
      const G4Track* secondary = particleChange->GetSecondary(i);
      
      // If neutral - eg, from photoelectric effect, just append
      if (secondary->GetDefinition()->GetPDGCharge() == 0) {
	survivors.push_back(new G4Track(*secondary)); 
      }
      else {
      // If charged, play Russian Roulette
	G4double random = G4UniformRand();
	
	if (random < inverseFactor) {
	  survivors.push_back(new G4Track(*secondary));
	}
      }
    }
    
    // Reconfigure particleChange to handle modified list of secondaries
    // The SetNumberOfSecondaries method deletes the old secondaries
    particleChange->SetNumberOfSecondaries(survivors.size());
    particleChange->SetSecondaryWeightByProcess(true);
    
    // New weight for charged secondaries
    G4double weight = track.GetWeight()*fFactor;
    
    // Add new secondaries to particle change 
    std::vector<G4Track*>::iterator iter = survivors.begin();
    
    while (iter != survivors.end()) {
      G4Track* myTrack = *iter;
      
      // Set weight to parent weight if secondary is a neutral
      if (myTrack->GetDefinition()->GetPDGCharge() == 0) weight = track.GetWeight();
      
      myTrack->SetWeight(weight);
      
      // particleChange takes ownership
      particleChange->AddSecondary(myTrack);
      
      iter++;
    }
    
  // Done
  return particleChange;
 
  }
}
#endif
