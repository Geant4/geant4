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
//  Gorad (Geant4 Open-source Radiation Analysis and Design)
//
//  Author : Makoto Asai (SLAC National Accelerator Laboratory)
//
//  Development of Gorad is funded by NASA Johnson Space Center (JSC)
//  under the contract NNJ15HK11B.
//
// ********************************************************************
//
// GRParallelWorldBiasingProcess.cc
//   A process that takes care of the geometry importance biasing.
//   This class assumes the existance of a parallel world dedicated
//   to the geometry importance biasing and biasing parameters are
//   associated to the region defined in that parallel world.
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#include "G4ios.hh"
#include "GRParallelWorldBiasingProcess.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Navigator.hh"
#include "G4PathFinder.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleChange.hh"
#include "G4ParticleChangeForNothing.hh"
#include "Randomize.hh"
#include "G4Region.hh"
#include "GRBiasingRegionInfo.hh"

GRParallelWorldBiasingProcess::
GRParallelWorldBiasingProcess(const G4String& processName,G4ProcessType theType)
:G4ParallelWorldProcess(processName,theType)
{
  particleChange = new G4ParticleChange();
  emptyParticleChange = new G4ParticleChangeForNothing();
}

GRParallelWorldBiasingProcess::~GRParallelWorldBiasingProcess()
{;}

G4VParticleChange* GRParallelWorldBiasingProcess::PostStepDoIt(
     const G4Track& track,
     const G4Step& step)
{ 
  fOldGhostTouchable = fGhostPostStepPoint->GetTouchableHandle();
  CopyStep(step);

  if(fOnBoundary)
  {
    fNewGhostTouchable = fPathFinder->CreateTouchableHandle(fNavigatorID);
  }
  else
  {
    fNewGhostTouchable = fOldGhostTouchable;
  }
    
  fGhostPreStepPoint->SetTouchableHandle(fOldGhostTouchable);
  fGhostPostStepPoint->SetTouchableHandle(fNewGhostTouchable);

  if(fNewGhostTouchable->GetVolume() == nullptr) 
  { 
    // world volume boundary
    emptyParticleChange->Initialize(track);
    return emptyParticleChange;
  }
  if(fNewGhostTouchable->GetVolume() == fOldGhostTouchable->GetVolume()) 
  {
    // stayed in the same volume
    emptyParticleChange->Initialize(track);
    return emptyParticleChange;
  }

  G4int preStepCopyNo = fOldGhostTouchable->GetReplicaNumber();
  G4int postStepCopyNo = fNewGhostTouchable->GetReplicaNumber();
  
  GRBiasingRegionInfo* biasInfo = static_cast<GRBiasingRegionInfo*>
    (fNewGhostTouchable->GetVolume()->GetLogicalVolume()->GetRegion()->GetUserInformation());
  auto probability = biasInfo->GetProbability();
  if(probability < 1.)
  {
    // skip biasing to avoid over biasing 
    G4double ran = G4UniformRand();
    if(ran > probability)
    {
      emptyParticleChange->Initialize(track);
      return emptyParticleChange;
    }
  }

  G4double trackWeight = track.GetWeight();
  particleChange->Initialize(track);

  auto biasFactor = biasInfo->GetBiasingFactor();
  if(biasFactor==1)
  {
    // not biasing
    emptyParticleChange->Initialize(track);
    return emptyParticleChange;
  }

  if(preStepCopyNo < postStepCopyNo)
  {
    // getting into more important volume, i.e. do splitting
    trackWeight /= biasFactor;
    particleChange->ProposeParentWeight(trackWeight);
    for(G4int iClone=1;iClone<biasFactor;iClone++)
    {
      G4Track* clonedTrack = new G4Track(track);
      clonedTrack->SetWeight(trackWeight);
      particleChange->AddSecondary(clonedTrack);
    }
    particleChange->SetSecondaryWeightByProcess(true);
  }
  else if(preStepCopyNo > postStepCopyNo)
  {
    // getting into less important volume, i.e. do Russian Rouletting
    G4double survivalRate = 1.0/biasFactor;
    G4double ran = G4UniformRand();
    if(ran > survivalRate)
    {
      // dead
      particleChange->ProposeTrackStatus(fStopAndKill);
    }
    else
    {
      // survive
      trackWeight *= biasFactor;
      particleChange->ProposeParentWeight(trackWeight);
    }
  } 
  else
  {
    // This should not happen!!
    G4ExceptionDescription ed;
    ed << "pre step point :  "<<fOldGhostTouchable->GetVolume()->GetName()<<" - copy no. "<<fOldGhostTouchable->GetReplicaNumber()<<"\n"
       << "post step point : "<<fNewGhostTouchable->GetVolume()->GetName()<<" - copy no. "<<fNewGhostTouchable->GetReplicaNumber();
    G4Exception("GRParallelWorldBiasingProcess::PostStepDoIt()","GoradProc0001",FatalException,ed);
  }

  return particleChange;
}

