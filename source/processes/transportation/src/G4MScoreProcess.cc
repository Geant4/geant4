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
// $Id: G4MScoreProcess.cc,v 1.8 2002-10-22 13:26:04 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
 : 
  G4VProcess(aName), 
  fScorer(aScorer),
  fKillTrack(false)
{
  G4VProcess::pParticleChange = new G4ParticleChange;
  if (!G4VProcess::pParticleChange) {
    G4std::G4Exception("ERROR:G4MScoreProcess::G4MScoreProcess new failed to create G4ParticleChange!");
  }
}

G4MScoreProcess::~G4MScoreProcess()
{
  delete pParticleChange;
}

G4double G4MScoreProcess::
PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				     G4double   previousStepSize,
				     G4ForceCondition* condition)
{
  *condition = Forced;
  return G4std::kInfinity;
}
  
G4VParticleChange * 
G4MScoreProcess::PostStepDoIt(const G4Track& aTrack, const G4Step &aStep)
{
  pParticleChange->Initialize(aTrack);

  if (aStep.GetStepLength() > G4std::kCarTolerance) {
    G4StepPoint *prepoint = aStep.GetPreStepPoint();
    G4StepPoint *postpoint = aStep.GetPostStepPoint();
  

    G4GeometryCell prekey(*(prepoint->GetPhysicalVolume()), 
			   prepoint->GetTouchable()->GetReplicaNumber());
    G4GeometryCell postkey(*(postpoint->GetPhysicalVolume()), 
			    postpoint->GetTouchable()->GetReplicaNumber());

    G4GeometryCellStep pstep(prekey, postkey);
    pstep.SetCrossBoundary(false);
    
    if (prekey != postkey) {
      pstep.SetCrossBoundary(true);
    } 
    fScorer.Score(aStep, pstep); 
  }

  if (fKillTrack) {
    fKillTrack = false;
    pParticleChange->SetStatusChange(fStopAndKill);
  }

  return G4VProcess::pParticleChange;
}

const G4String &G4MScoreProcess::GetName() const {
  return theProcessName;
}


G4double G4MScoreProcess::
AlongStepGetPhysicalInteractionLength(const G4Track&,
				      G4double  ,
				      G4double  ,
				      G4double& ,
				      G4GPILSelection*) {
  return -1.0;
}

G4double G4MScoreProcess::
AtRestGetPhysicalInteractionLength(const G4Track&,
				   G4ForceCondition*) {
  return -1.0;
}

G4VParticleChange* G4MScoreProcess::AtRestDoIt(const G4Track&,
					       const G4Step&) {
  return 0;
}

G4VParticleChange* G4MScoreProcess::AlongStepDoIt(const G4Track&,
						  const G4Step&) {
  return 0;
}

void G4MScoreProcess::KillTrack() const{
  fKillTrack = true;
}
