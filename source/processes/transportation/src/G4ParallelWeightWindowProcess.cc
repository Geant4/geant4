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
// $Id: G4ParallelWeightWindowProcess.cc,v 1.10 2003/11/26 14:51:50 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ParallelWeightWindowProcess.cc
//
// ----------------------------------------------------------------------

#include "G4Types.hh"
#include <strstream>
#include "G4ParallelWeightWindowProcess.hh"
#include "G4VWeightWindowExaminer.hh"
#include "G4VTrackTerminator.hh"
#include "G4SamplingPostStepAction.hh"
#include "G4VParallelStepper.hh"

G4ParallelWeightWindowProcess::G4ParallelWeightWindowProcess(
                         const G4VWeightWindowExaminer &aWeightWindowExaminer,
                               G4VParallelStepper &aStepper,
                         const G4VTrackTerminator *TrackTerminator,
                               G4PlaceOfAction placeOfAction,
                         const G4String &aName) 
  : G4VProcess(aName),
    fParticleChange(new G4ParticleChange),
    fWeightWindowExaminer(aWeightWindowExaminer),
    fStepper(aStepper),
    fPostStepAction(0),
    fPlaceOfAction(placeOfAction)
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

G4ParallelWeightWindowProcess::~G4ParallelWeightWindowProcess()
{
  delete fParticleChange;
  delete fPostStepAction;
}

G4double G4ParallelWeightWindowProcess::
PostStepGetPhysicalInteractionLength(const G4Track& ,
                                     G4double  ,
                                     G4ForceCondition* condition)
{
  *condition = Forced;
  return kInfinity;
}

G4VParticleChange *
G4ParallelWeightWindowProcess::PostStepDoIt(const G4Track& aTrack,
                                            const G4Step &aStep)
{
  if (aTrack.GetTrackStatus()==fStopAndKill)
  {
    G4cout << "WARNING - G4ParallelWeightWindowProcess::PostStepDoIt()"
           << "          StopAndKill track !" << G4endl;
  }
  
  fParticleChange->Initialize(aTrack);
  
  if (aStep.GetPostStepPoint()->GetStepStatus() != fGeomBoundary)
  {
    if (!( fStepper.GetPStep().GetCrossBoundary()
        && (fPlaceOfAction == onCollision) ) )
    {
      // get new weight and number of clones
      G4Nsplit_Weight nw(fWeightWindowExaminer.
                         Examine(aTrack.GetWeight(),
                                 aTrack.
                                 GetKineticEnergy()));
      fPostStepAction->DoIt(aTrack, fParticleChange, nw);
    }
  }
  return fParticleChange;
}
  
void G4ParallelWeightWindowProcess::Error(const G4String &m)
{
  G4Exception("G4ParallelWeightWindowProcess::Error()",
              "ProgramError", FatalException, m);
}

void G4ParallelWeightWindowProcess::KillTrack() const
{
  fParticleChange->SetStatusChange(fStopAndKill);
}

const G4String &G4ParallelWeightWindowProcess::GetName() const
{
  return theProcessName;
}

G4double G4ParallelWeightWindowProcess::
AlongStepGetPhysicalInteractionLength(const G4Track&,
                                      G4double  ,
                                      G4double  ,
                                      G4double& ,
                                      G4GPILSelection*)
{
  return -1.0;
}

G4double G4ParallelWeightWindowProcess::
AtRestGetPhysicalInteractionLength(const G4Track& ,
                                   G4ForceCondition*) 
{
  return -1.0;
}
  
G4VParticleChange* G4ParallelWeightWindowProcess::
AtRestDoIt(const G4Track&, const G4Step&) 
{
  return 0;
}

G4VParticleChange* G4ParallelWeightWindowProcess::
AlongStepDoIt(const G4Track&, const G4Step&)
{
  return 0;
}
