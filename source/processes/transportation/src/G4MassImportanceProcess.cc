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
// $Id: G4MassImportanceProcess.cc,v 1.17 2004/10/19 00:59:39 kurasige Exp $
// GEANT4 tag $Name: geant4-07-01 $
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

