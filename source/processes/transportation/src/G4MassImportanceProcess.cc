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
// $Id: G4MassImportanceProcess.cc,v 1.9 2002-10-16 16:27:00 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

G4MassImportanceProcess::
G4MassImportanceProcess(const G4VImportanceAlgorithm &aImportanceAlgorithm,
			const G4VIStore &aIstore,
			const G4VTrackTerminator *TrackTerminator,
			const G4String &aName)
 : G4VProcess(aName),
   fParticleChange(new G4ParticleChange),
   fTrackTerminator(TrackTerminator ? TrackTerminator : this),
   fImportanceAlgorithm(aImportanceAlgorithm),
   fImportanceFinder(aIstore),
   fImportancePostStepDoIt(*fTrackTerminator)
{
  if (!fParticleChange) {
    G4std::G4Exception("ERROR:G4MassImportanceProcess::G4MassImportanceProcess: new failed to create G4ParticleChange!");
  }
  G4VProcess::pParticleChange = fParticleChange;
}

G4MassImportanceProcess::~G4MassImportanceProcess()
{
  delete fParticleChange;
}

G4double G4MassImportanceProcess::
PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				     G4double   previousStepSize,
				     G4ForceCondition* condition)
{
  *condition = Forced;
  return G4std::kInfinity;
}
  
G4VParticleChange *
G4MassImportanceProcess::PostStepDoIt(const G4Track &aTrack,
				      const G4Step &aStep)
{
  fParticleChange->Initialize(aTrack);
  if (aStep.GetPostStepPoint()->GetStepStatus() == fGeomBoundary
      && aStep.GetStepLength() > G4std::kCarTolerance) {
    if (aTrack.GetTrackStatus()==fStopAndKill) {
      G4std::G4cout << "G4MassImportanceProcess::PostStepDoIt StopAndKill" << G4endl;
    }

    G4StepPoint *prepoint = aStep.GetPreStepPoint();
    G4StepPoint *postpoint = aStep.GetPostStepPoint();
  

    G4GeometryCell prekey(*(prepoint->GetPhysicalVolume()), 
			 prepoint->GetTouchable()->GetReplicaNumber());
    G4GeometryCell postkey(*(postpoint->GetPhysicalVolume()), 
			  postpoint->GetTouchable()->GetReplicaNumber());

    G4Nsplit_Weight nw = fImportanceAlgorithm.
      Calculate(fImportanceFinder.GetImportance(prekey),
		fImportanceFinder.GetImportance(postkey), 
		aTrack.GetWeight());
    fImportancePostStepDoIt.DoIt(aTrack, fParticleChange, nw);
  }
  return fParticleChange;
}

void G4MassImportanceProcess::KillTrack() const {
  fParticleChange->SetStatusChange(fStopAndKill);
}

const G4String &G4MassImportanceProcess::GetName() const {
  return theProcessName;
}

G4double G4MassImportanceProcess::
AlongStepGetPhysicalInteractionLength(const G4Track&,
				      G4double  ,
				      G4double  ,
				      G4double& ,
				      G4GPILSelection*) {
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
AlongStepDoIt(const G4Track&, const G4Step&) {
  return 0;
}

