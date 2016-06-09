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
// $Id: G4ParallelImportanceProcess.cc,v 1.18 2005/11/21 21:41:29 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ParallelImportanceProcess.cc
//
// ----------------------------------------------------------------------

#include "G4Types.hh"
#include "G4ParallelImportanceProcess.hh"
#include "G4VImportanceSplitExaminer.hh"
#include "G4VTrackTerminator.hh"
#include "G4SamplingPostStepAction.hh"

G4ParallelImportanceProcess::G4ParallelImportanceProcess(
                 const G4VImportanceSplitExaminer &aImportanceSplitExaminer,
                       G4VPGeoDriver &pgeodriver,
                       G4VParallelStepper &aStepper, 
                 const G4VTrackTerminator *TrackTerminator,
                 const G4String &aName)
 : G4ParallelTransport(pgeodriver, aStepper, aName),
   fParticleChange(G4ParallelTransport::fParticleChange),
   fImportanceSplitExaminer(aImportanceSplitExaminer),
   fPostStepAction(0)
{
  if (TrackTerminator)
  {
    fPostStepAction = new G4SamplingPostStepAction(*TrackTerminator);
  }
  else
  {
    fPostStepAction = new G4SamplingPostStepAction(*this);
  }

}

G4ParallelImportanceProcess::~G4ParallelImportanceProcess()
{
  delete fPostStepAction;
}

G4VParticleChange *G4ParallelImportanceProcess::
PostStepDoIt(const G4Track& aTrack, const G4Step &aStep)
{
  if (aTrack.GetTrackStatus()==fStopAndKill)
  {
    G4cout << "WARNING - G4ParallelImportanceProcess::PostStepDoIt()"
           << "          StopAndKill track !" << G4endl;
  }
  G4ParallelTransport::PostStepDoIt(aTrack, aStep);

  // get new weight and number of clones
  G4Nsplit_Weight nw(fImportanceSplitExaminer.Examine(aTrack.GetWeight()));

  fPostStepAction->DoIt(aTrack, fParticleChange, nw);
  return fParticleChange;
}
  
void G4ParallelImportanceProcess::Error(const G4String &m)
{
  G4Exception("G4ParallelImportanceProcess::Error()",
              "ProgramError", FatalException, m);
}

void G4ParallelImportanceProcess::KillTrack() const
{
  fParticleChange->ProposeTrackStatus(fStopAndKill);
}

const G4String &G4ParallelImportanceProcess::GetName() const
{
  return G4ParallelTransport::GetProcessName();
}
