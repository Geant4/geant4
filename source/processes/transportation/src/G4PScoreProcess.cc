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
//
// $Id: G4PScoreProcess.cc,v 1.3 2002-08-13 10:07:47 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4PScoreProcess.cc
//
// ----------------------------------------------------------------------

#include "G4PScoreProcess.hh"
#include "G4VPScorer.hh"
#include "G4PStep.hh"
#include "G4VParallelStepper.hh"

G4PScoreProcess::G4PScoreProcess(G4VParallelStepper &astepper,
				 G4VPScorer &aScorer,
				 const G4String &aName)
 : 
  G4VProcess(aName), 
  fPstepper(astepper),
  fScorer(aScorer),
  fKillTrack(false)
{
  G4VProcess::pParticleChange = new G4ParticleChange;
}

G4PScoreProcess::~G4PScoreProcess()
{
  delete pParticleChange;
}

G4double G4PScoreProcess::
PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				     G4double   previousStepSize,
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
  if (fKillTrack) {
    fKillTrack = false;
    pParticleChange->SetStatusChange(fStopAndKill);
  }
  return G4VProcess::pParticleChange;
}

