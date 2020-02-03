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
/// \file GB06BOptnSplitAndKillByImportance.cc
/// \brief Implementation of the GB06BOptnSplitAndKillByImportance class

#include "GB06BOptnSplitAndKillByImportance.hh"
#include "Randomize.hh"


#include "G4BiasingProcessInterface.hh"
#include "G4ParallelGeometriesLimiterProcess.hh"
#include "G4BiasingProcessSharedData.hh"
#include "G4ProcessManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB06BOptnSplitAndKillByImportance::GB06BOptnSplitAndKillByImportance(G4String name)
: G4VBiasingOperation (name),
  fParallelWorldIndex ( -1 ),
  fBiasingSharedData  ( nullptr ),
  fParticleChange(),
  fDummyParticleChange()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB06BOptnSplitAndKillByImportance::~GB06BOptnSplitAndKillByImportance()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double GB06BOptnSplitAndKillByImportance::
DistanceToApplyOperation( const G4Track*,
                          G4double,
                          G4ForceCondition* condition)
{

  // -- Remember the touchable history (ie geometry state) at the beginning of the step:
  // -- Start by getting the process handling the step limitation in parallel
  // -- geometries for the generic biasing:
  auto biasingLimiterProcess =  fBiasingSharedData->GetParallelGeometriesLimiterProcess();
  fPreStepTouchableHistory   =
    biasingLimiterProcess
    ->GetNavigator( fParallelWorldIndex )  // -- get the navigator of the geometry...
    ->CreateTouchableHistoryHandle();      // -- ...and collect the geometry state.
  
  // -- We request to be "forced" : meaning GenerateBiasingFinalState(...) will be called
  // -- anyway at the PostStepDoIt(...) stage (ie, when the track will have been moved to
  // -- its end point position) and, there, we will have to handle the decision to
  // -- split/kill if we are on a volume boundary, or do nothing, if we're not:
  *condition = Forced;

  // -- As we're forced, we make this returned length void:
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* 
GB06BOptnSplitAndKillByImportance::GenerateBiasingFinalState( const G4Track* track,
                                                              const G4Step*       )
{
  // -- Given we used the "Forced" condition, this method is always called.
  // -- We check the status of the step in the parallel geometry, and apply
  // -- splitting/killing if the step has been limited by the geometry.

  // -- We first check if we limit the step in the parallel geometry:
  G4bool isLimiting = fBiasingSharedData
    ->GetParallelGeometriesLimiterProcess()
    ->GetIsLimiting( fParallelWorldIndex );

  // -- if not limiting, we leave the track unchanged:
  if ( !isLimiting )
    {
      fDummyParticleChange.Initialize( *track );
      return &fDummyParticleChange;
    }
  
  // -- We collect the geometry state at the end point step:
  // -- Note this is the same call than in the DistanceToApplyOperation(...) for the
  // -- fPreStepTouchableHistory, but the navigator has pushed the track in the next
  // -- volume since then (even if the track is still on the boundary), and hence the
  // -- geometry state has changed.
  auto biasingLimiterProcess = fBiasingSharedData->GetParallelGeometriesLimiterProcess();
  fPostStepTouchableHistory  =
    biasingLimiterProcess
    ->GetNavigator( fParallelWorldIndex )
    ->CreateTouchableHistoryHandle();

  
  // -- We verify we are still in the "same" physical volume:
  // -- remember : using a replica, we have "one" physical volume
  // -- but representing many different placements, distinguished
  // -- by replica number. By checking we are in the same physical
  // -- volume, we verify the end step point has not reached the
  // -- world volume of the parallel geometry:
  if ( fPreStepTouchableHistory ->GetVolume() !=
       fPostStepTouchableHistory->GetVolume()    )
    {
      // -- the track is leaving the volumes in which biasing is applied,
      // -- we leave this track unchanged:
      fDummyParticleChange.Initialize( *track );
      return &fDummyParticleChange;
    }

  // -------------------------------------------------------------------------------------
  // -- At this stage, we know we have a particle crossing a boundary between two slices,
  // -- each of this slice has a well defined importance : we apply the biasing.
  // -- We will split if the track is entering a volume with higher importance, and
  // -- kill (applying Russian roulette) in the other case.
  // -------------------------------------------------------------------------------------
  
  // -- We start by getting the replica numbers: 
  G4int  preReplicaNumber =  fPreStepTouchableHistory->GetReplicaNumber();
  G4int postReplicaNumber = fPostStepTouchableHistory->GetReplicaNumber();

  // -- and get the related importance we defined in the importance map:
  G4int  preImportance = (*fImportanceMap)[ preReplicaNumber];
  G4int postImportance = (*fImportanceMap)[postReplicaNumber];


  // -- Get track weight:
  G4double initialWeight = track->GetWeight();
  
  // -- Initialize the "particle change" (it will communicate the new track state to
  // -- the tracking):
  fParticleChange.Initialize(*track);

  if ( postImportance > preImportance )
    {
      // -- We split :
      G4int splittingFactor = postImportance/preImportance;
      
      // Define the tracks weight:
      G4double weightOfTrack = initialWeight/splittingFactor;
      
      // Ask currect track weight to be changed to the new value:
      fParticleChange.ProposeParentWeight( weightOfTrack );
      // Now we clone this track (this is the actual splitting):
      // we will then have the primary and clone of it, hence the
      // splitting by a factor 2:
      G4Track* clone = new G4Track( *track );
      clone->SetWeight( weightOfTrack );
      fParticleChange.AddSecondary( clone );
      // -- Below's call added for safety & illustration : inform particle change to not
      // -- modify the clone (ie : daughter) weight to make it that of the primary.
      // -- Here this call is not mandatory as both tracks have same weights.
      fParticleChange.SetSecondaryWeightByProcess(true);
    }
  else
    {
      // -- We apply Russian roulette:
      G4double survivingProbability = G4double(postImportance) / G4double(preImportance);
      
      // Shoot a random number (in ]0,1[ segment):
      G4double random = G4UniformRand();
      if ( random > survivingProbability )
        {
          // We ask for the the track to be killed:
          fParticleChange.ProposeTrackStatus(fStopAndKill);
        }
      else
        {
          // In this case, the track survives. We change its weight
          // to conserve weight among killed and survival tracks:
          fParticleChange.ProposeParentWeight( initialWeight/survivingProbability );
        }
    }

  return &fParticleChange;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
