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
// $Id: G4MassWeightWindowProcess.cc,v 1.4 2003/11/26 14:51:49 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4MassWeightWindowProcess.cc
//
// ----------------------------------------------------------------------

#include "G4MassWeightWindowProcess.hh"
#include "G4VWeightWindowAlgorithm.hh"
#include "G4GeometryCell.hh"
#include "G4SamplingPostStepAction.hh"
#include "G4VTrackTerminator.hh"
#include "G4PlaceOfAction.hh"
#include "G4VWeightWindowStore.hh"

G4MassWeightWindowProcess::G4MassWeightWindowProcess(
                    const G4VWeightWindowAlgorithm &aWeightWindowAlgorithm,
                    const G4VWeightWindowStore &aWWStore,
                    const G4VTrackTerminator *TrackTerminator,
                          G4PlaceOfAction placeOfAction,
                    const G4String &aName)
 : G4VProcess(aName),
   fParticleChange(new G4ParticleChange),
   fWeightWindowAlgorithm(aWeightWindowAlgorithm),
   fWeightWindowStore(aWWStore),
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
  if (!fParticleChange)
  {
    G4Exception("G4MassWeightWindowProcess::G4MassWeightWindowProcess()",
                "FatalError", FatalException,
                "Failed allocation of G4ParticleChange !");
  }
  G4VProcess::pParticleChange = fParticleChange;
}

G4MassWeightWindowProcess::~G4MassWeightWindowProcess()
{
  delete fPostStepAction;
  delete fParticleChange;
}

G4double G4MassWeightWindowProcess::
PostStepGetPhysicalInteractionLength(const G4Track& ,
                                     G4double   ,
                                     G4ForceCondition* condition)
{
  *condition = Forced;
  return kInfinity;
}
  
G4VParticleChange *
G4MassWeightWindowProcess::PostStepDoIt(const G4Track &aTrack,
                                        const G4Step &aStep)
{
  fParticleChange->Initialize(aTrack);
  if (aStep.GetStepLength() > kCarTolerance)
  {
    if ( ( fPlaceOfAction == onBoundaryAndCollision)
      || ( (fPlaceOfAction == onBoundary) && 
           (aStep.GetPostStepPoint()->GetStepStatus() == fGeomBoundary) )
      || ( (fPlaceOfAction == onCollision) && 
           (aStep.GetPostStepPoint()->GetStepStatus() != fGeomBoundary) ) )
    {

      G4StepPoint *postpoint = aStep.GetPostStepPoint();

      G4GeometryCell postCell(*(postpoint->GetPhysicalVolume()), 
                             postpoint->GetTouchable()->GetReplicaNumber());

      G4Nsplit_Weight nw =
        fWeightWindowAlgorithm.Calculate(aTrack.GetWeight(),
                  fWeightWindowStore.GetLowerWeitgh(postCell,
                                                    aTrack.GetKineticEnergy()));
      fPostStepAction->DoIt(aTrack, fParticleChange, nw);
    }
  }
  return fParticleChange;
}

void G4MassWeightWindowProcess::KillTrack() const
{
  fParticleChange->SetStatusChange(fStopAndKill);
}

const G4String &G4MassWeightWindowProcess::GetName() const
{
  return G4VProcess::GetProcessName();
}

G4double G4MassWeightWindowProcess::
AlongStepGetPhysicalInteractionLength(const G4Track&,
                                      G4double  ,
                                      G4double  ,
                                      G4double& ,
                                      G4GPILSelection*)
{
  return -1.0;
}

G4double G4MassWeightWindowProcess::
AtRestGetPhysicalInteractionLength(const G4Track& ,
                                   G4ForceCondition*) 
{
  return -1.0;
}
  
G4VParticleChange* G4MassWeightWindowProcess::
AtRestDoIt(const G4Track&, const G4Step&) 
{
  return 0;
}

G4VParticleChange* G4MassWeightWindowProcess::
AlongStepDoIt(const G4Track&, const G4Step&)
{
  return 0;
}
