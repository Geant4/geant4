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
// $Id: BremSplittingProcess.cc,v 1.3 2007-06-21 15:04:20 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Jane Tinslay, May 2007
//
#include "BremSplittingProcess.hh"

#include "ConfigData.hh"
#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include <assert.h>
#include <vector>

BremSplittingProcess::BremSplittingProcess() {}

BremSplittingProcess::~BremSplittingProcess() {}

G4VParticleChange* 
BremSplittingProcess::PostStepDoIt(const G4Track& track, const G4Step& step)
{
  // Just do regular processing if brem splitting is not activated
  G4VParticleChange* particleChange(0);
  unsigned nSplit = ConfigData::BremSplitting::GetFactor();

  if (!ConfigData::BremSplitting::GetActivation()) {
    particleChange = pRegProcess->PostStepDoIt(track, step);
    assert (0 != particleChange);

    ConfigData::BremSplitting::AddSecondaries(particleChange->GetNumberOfSecondaries());

    return particleChange;
  }
  
  // Do brem splitting
  assert (nSplit > 0);

  unsigned i(0);
  G4double weight = track.GetWeight()/nSplit;
  
  // Secondary store
  std::vector<G4Track*> secondaries;
  secondaries.reserve(nSplit);
    
  // Loop over PostStepDoIt method to generate multiple secondaries.
  for (i=0; i<nSplit; i++) {    
    particleChange = pRegProcess->PostStepDoIt(track, step);

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

  ConfigData::BremSplitting::AddSecondaries(secondaries.size());

  return particleChange;
}
