//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: BremSplittingProcess.cc,v 1.2 2007-05-24 22:05:50 tinslay Exp $
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
