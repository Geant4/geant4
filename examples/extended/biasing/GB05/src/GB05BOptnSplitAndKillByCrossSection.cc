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
/// \file GB05BOptnSplitAndKillByCrossSection.cc
/// \brief Implementation of the GB05BOptnSplitAndKillByCrossSection class

#include "Randomize.hh"
#include "GB05BOptnSplitAndKillByCrossSection.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB05BOptnSplitAndKillByCrossSection::GB05BOptnSplitAndKillByCrossSection(G4String name)
: G4VBiasingOperation(name),
  fParticleChange(),
  fInteractionLength(-1.0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB05BOptnSplitAndKillByCrossSection::~GB05BOptnSplitAndKillByCrossSection()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double GB05BOptnSplitAndKillByCrossSection::
DistanceToApplyOperation( const G4Track*,
                          G4double,
                          G4ForceCondition* condition)
{
  *condition = NotForced;

  // -- Sample the exponential law using the total interaction length of processes
  // -- to couterbalance for:
  G4double proposedStepLength =  -std::log( G4UniformRand() ) * fInteractionLength;
  
  return proposedStepLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* 
GB05BOptnSplitAndKillByCrossSection::GenerateBiasingFinalState( const G4Track* track,
                                                                const G4Step*       )
{
  
  // -- This method is called if we have limited the step.
  // -- We hence make the splitting or killing.

  // Get track weight:
  G4double initialWeight = track->GetWeight();
  
  // The "particle change" is the object to be used to communicate to
  // the tracking the update of the primary state and/or creation
  // secondary tracks.
  fParticleChange.Initialize(*track);

  // -- Splitting and killing factors.
  // -- They are taken the same, but the killing factor can be make bigger.
  G4double splittingFactor =  2.0;
  G4double   killingFactor =  2.0;


  if ( track->GetMomentumDirection().z() > 0 )
    {
      // -- We split if the track is moving forward:
      
      // Define the tracks weight:
      G4double weightOfTrack = initialWeight/splittingFactor;
      
      // Ask currect track weight to be changed to new value:
      fParticleChange.ProposeParentWeight( weightOfTrack );
      // Now make clones of this track (this is the actual splitting):
      // we will then have the primary and clone of it, hence the
      // splitting by a factor 2:
      G4Track* clone = new G4Track( *track );
      clone->SetWeight( weightOfTrack );
      fParticleChange.AddSecondary( clone );
      // -- Below's call added for safety & illustration : inform particle change to not
      // -- modify the clone (ie : daughter) weight to male it that of the
      // -- primary. Here call not mandatory and both tracks have same weights.
      fParticleChange.SetSecondaryWeightByProcess(true);
    }
  else
    {
      // -- We apply Russian roulette if the track is moving backward:
      
      // Shoot a random number (in ]0,1[ segment):
      G4double random = G4UniformRand();
      G4double killingProbability = 1.0 - 1.0/killingFactor;
      if ( random < killingProbability )
        {
          // We ask for the the track to be killed:
          fParticleChange.ProposeTrackStatus(fStopAndKill);
        }
      else
        {
          // In this case, the track survives. We change its weight
          // to conserve weight among killed and survival tracks:
          fParticleChange.ProposeParentWeight( initialWeight*killingFactor );
        }
    }
  
  return &fParticleChange;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
