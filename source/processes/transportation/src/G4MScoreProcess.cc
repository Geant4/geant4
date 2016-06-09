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
// $Id: G4MScoreProcess.cc,v 1.14 2006/11/13 16:19:14 japost Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4MScoreProcess.cc
//
// ----------------------------------------------------------------------

#include "G4MScoreProcess.hh"
#include "G4VScorer.hh"
#include "G4GeometryCellStep.hh"
#include "G4VParallelStepper.hh"

G4MScoreProcess::G4MScoreProcess(G4VScorer &aScorer,
                                 const G4String &aName)
 : G4VProcess(aName), 
   fScorer(aScorer),
   fKillTrack(false)
{
  G4VProcess::pParticleChange = new G4ParticleChange;
  if (!G4VProcess::pParticleChange)
  {
    G4Exception("G4MScoreProcess::G4MScoreProcess()", "Fatal Error",
                FatalException, "Failed to allocate G4ParticleChange !");
  }
}

G4MScoreProcess::~G4MScoreProcess()
{
  delete pParticleChange;
}

G4double G4MScoreProcess::
PostStepGetPhysicalInteractionLength(const G4Track&,
                                     G4double  ,
                                     G4ForceCondition* condition)
{
  *condition = Forced;
  return kInfinity;
}
  
G4VParticleChange * 
G4MScoreProcess::PostStepDoIt(const G4Track& aTrack, const G4Step &aStep)
{
  pParticleChange->Initialize(aTrack);

  if (aStep.GetStepLength() > kCarTolerance)
  {
    G4StepPoint *prepoint = aStep.GetPreStepPoint();
    G4StepPoint *postpoint = aStep.GetPostStepPoint();

    G4GeometryCell prekey(*(prepoint->GetPhysicalVolume()), 
                           prepoint->GetTouchable()->GetReplicaNumber());
    G4GeometryCell postkey(*(postpoint->GetPhysicalVolume()), 
                            postpoint->GetTouchable()->GetReplicaNumber());

    G4GeometryCellStep pstep(prekey, postkey);
    pstep.SetCrossBoundary(false);
    
    if (prekey != postkey)
    {
      pstep.SetCrossBoundary(true);
    } 
    fScorer.Score(aStep, pstep); 
  }

  if (fKillTrack)
  {
    fKillTrack = false;
    pParticleChange->ProposeTrackStatus(fStopAndKill);
  }

  return G4VProcess::pParticleChange;
}

const G4String &G4MScoreProcess::GetName() const
{
  return theProcessName;
}


G4double G4MScoreProcess::
AlongStepGetPhysicalInteractionLength(const G4Track&,
                                      G4double  ,
                                      G4double  ,
                                      G4double& ,
                                      G4GPILSelection*)
{
  return -1.0;
}

G4double G4MScoreProcess::
AtRestGetPhysicalInteractionLength(const G4Track&,
                                   G4ForceCondition*)
{
  return -1.0;
}

G4VParticleChange* G4MScoreProcess::AtRestDoIt(const G4Track&,
                                               const G4Step&)
{
  return 0;
}

G4VParticleChange* G4MScoreProcess::AlongStepDoIt(const G4Track&,
                                                  const G4Step&)
{
  return 0;
}

void G4MScoreProcess::KillTrack() const
{
  fKillTrack = true;
}
