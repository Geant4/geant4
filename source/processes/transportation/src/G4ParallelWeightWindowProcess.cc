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
// $Id: G4ParallelWeightWindowProcess.cc,v 1.13 2006/06/29 21:12:20 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ParallelWeightWindowProcess.cc
//
// ----------------------------------------------------------------------

#include "G4Types.hh"
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
  fParticleChange->ProposeTrackStatus(fStopAndKill);
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
