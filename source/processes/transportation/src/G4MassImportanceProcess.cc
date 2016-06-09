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
// $Id: G4MassImportanceProcess.cc,v 1.17.2.1 2006/06/29 21:11:02 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4MassImportanceProcess.cc
//
// ----------------------------------------------------------------------

#include "G4MassImportanceProcess.hh"
#include "G4VImportanceAlgorithm.hh"
#include "G4GeometryCell.hh"
#include "G4SamplingPostStepAction.hh"
#include "G4VTrackTerminator.hh"
#include "G4VIStore.hh"

G4MassImportanceProcess::
G4MassImportanceProcess(const G4VImportanceAlgorithm &aImportanceAlgorithm,
                        const G4VIStore &aIstore,
                        const G4VTrackTerminator *TrackTerminator,
                        const G4String &aName)
 : G4VProcess(aName),
   fParticleChange(new G4ParticleChange),
   fImportanceAlgorithm(aImportanceAlgorithm),
   fIStore(aIstore),
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
  if (!fParticleChange)
  {
    G4Exception("G4MassImportanceProcess::G4MassImportanceProcess()",
                "FatalError", FatalException,
                "Failed allocation of G4ParticleChange !");
  }
  G4VProcess::pParticleChange = fParticleChange;
}

G4MassImportanceProcess::~G4MassImportanceProcess()
{
  delete fPostStepAction;
  delete fParticleChange;
}

G4double G4MassImportanceProcess::
PostStepGetPhysicalInteractionLength(const G4Track& ,
                                     G4double   ,
                                     G4ForceCondition* condition)
{
  *condition = Forced;
  return kInfinity;
}
  
G4VParticleChange *
G4MassImportanceProcess::PostStepDoIt(const G4Track &aTrack,
                                      const G4Step &aStep)
{
  fParticleChange->Initialize(aTrack);
  if ( (aStep.GetPostStepPoint()->GetStepStatus() == fGeomBoundary)
    && (aStep.GetStepLength() > kCarTolerance) )
  {
    if (aTrack.GetTrackStatus()==fStopAndKill)
    {
      G4cout << "WARNING - G4MassImportanceProcess::PostStepDoIt()"
             << "          StopAndKill track." << G4endl;
    }

    G4StepPoint *prepoint  = aStep.GetPreStepPoint();
    G4StepPoint *postpoint = aStep.GetPostStepPoint();

    G4GeometryCell prekey(*(prepoint->GetPhysicalVolume()), 
                         prepoint->GetTouchable()->GetReplicaNumber());
    G4GeometryCell postkey(*(postpoint->GetPhysicalVolume()), 
                          postpoint->GetTouchable()->GetReplicaNumber());

    G4Nsplit_Weight nw = fImportanceAlgorithm.
      Calculate(fIStore.GetImportance(prekey),
                fIStore.GetImportance(postkey), 
                aTrack.GetWeight());
    fPostStepAction->DoIt(aTrack, fParticleChange, nw);
  }
  return fParticleChange;
}

void G4MassImportanceProcess::KillTrack() const
{
  fParticleChange->ProposeTrackStatus(fStopAndKill);
}

const G4String &G4MassImportanceProcess::GetName() const
{
  return theProcessName;
}

G4double G4MassImportanceProcess::
AlongStepGetPhysicalInteractionLength(const G4Track&,
                                      G4double  ,
                                      G4double  ,
                                      G4double& ,
                                      G4GPILSelection*)
{
  return -1.0;
}

G4double G4MassImportanceProcess::
AtRestGetPhysicalInteractionLength(const G4Track& ,
                                   G4ForceCondition*) 
{
  return -1.0;
}
  
G4VParticleChange* G4MassImportanceProcess::
AtRestDoIt(const G4Track&, const G4Step&) 
{
  return 0;
}

G4VParticleChange* G4MassImportanceProcess::
AlongStepDoIt(const G4Track&, const G4Step&)
{
  return 0;
}

