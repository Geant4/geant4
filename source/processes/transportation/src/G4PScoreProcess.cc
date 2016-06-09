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
// $Id: G4PScoreProcess.cc,v 1.9.2.1 2006/06/29 21:12:04 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4PScoreProcess.cc
//
// ----------------------------------------------------------------------

#include "G4PScoreProcess.hh"
#include "G4VScorer.hh"
#include "G4GeometryCellStep.hh"
#include "G4VParallelStepper.hh"

G4PScoreProcess::G4PScoreProcess(G4VParallelStepper &astepper,
                                 G4VScorer &aScorer,
                                 const G4String &aName)
 : G4VProcess(aName), 
   fPstepper(astepper),
   fScorer(aScorer),
   fKillTrack(false)
{
  G4VProcess::pParticleChange = new G4ParticleChange;
  if (!G4VProcess::pParticleChange)
  {
    G4Exception("G4PScoreProcess::G4PScoreProcess()", "FatalError",
                FatalException, "Failed to allocate G4ParticleChange !");
  }
}

G4PScoreProcess::~G4PScoreProcess()
{
  delete pParticleChange;
}

G4double G4PScoreProcess::
PostStepGetPhysicalInteractionLength(const G4Track &,
                                     G4double ,
                                     G4ForceCondition* condition)
{
  *condition = Forced;
  return kInfinity;
}
  
G4VParticleChange * 
G4PScoreProcess::PostStepDoIt(const G4Track& aTrack, const G4Step &aStep)
{
  pParticleChange->Initialize(aTrack);
  fScorer.Score(aStep, fPstepper.GetPStep()); 
  if (fKillTrack)
  {
    fKillTrack = false;
    pParticleChange->ProposeTrackStatus(fStopAndKill);
  }
  return G4VProcess::pParticleChange;
}

void G4PScoreProcess::KillTrack() const
{
  fKillTrack = true;
}

const G4String &G4PScoreProcess::GetName() const
{
  return theProcessName;
}

G4double G4PScoreProcess::
AlongStepGetPhysicalInteractionLength(const G4Track&,
                                      G4double  ,
                                      G4double  ,
                                      G4double& ,
                                      G4GPILSelection*)
{
  return -1.0;
}

G4double G4PScoreProcess::
AtRestGetPhysicalInteractionLength(const G4Track& , G4ForceCondition*)
{
  return -1.0;
}

G4VParticleChange* G4PScoreProcess::
AtRestDoIt(const G4Track&, const G4Step&)
{
  return 0;
}

G4VParticleChange* G4PScoreProcess::
AlongStepDoIt(const G4Track&, const G4Step&)
{
  return 0;
}
