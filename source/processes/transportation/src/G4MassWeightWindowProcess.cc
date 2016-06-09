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
// $Id: G4MassWeightWindowProcess.cc,v 1.5.2.1 2006/06/29 21:11:04 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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
  fParticleChange->ProposeTrackStatus(fStopAndKill);
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
