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
// $Id: G4ParallelWWnTransportProcess.cc,v 1.1 2003-08-19 15:18:22 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ParallelWWnTransportProcess.cc
//
// ----------------------------------------------------------------------

#include "G4Types.hh"
#include <strstream>
#include "G4ParallelWWnTransportProcess.hh"
#include "G4VWeightWindowExaminer.hh"
#include "G4VTrackTerminator.hh"
#include "G4ImportancePostStepDoIt.hh"

G4ParallelWWnTransportProcess::
G4ParallelWWnTransportProcess(const G4VWeightWindowExaminer 
			      &aWeightWindowExaminer,
			      G4VPGeoDriver &pgeodriver, 
			      G4VParallelStepper &aStepper,
			      const G4VTrackTerminator *TrackTerminator,
			      const G4String &aName) 
  : 
  G4ParallelTransport(pgeodriver, aStepper, aName),
  fParticleChange(G4ParallelTransport::fParticleChange),
  fWeightWindowExaminer(aWeightWindowExaminer),
  fImportancePostStepDoIt(0)
{
  if (TrackTerminator) {
    fImportancePostStepDoIt = new G4ImportancePostStepDoIt(*TrackTerminator);
  }
  else {
    fImportancePostStepDoIt = new G4ImportancePostStepDoIt(*this);
  }

}

G4ParallelWWnTransportProcess::~G4ParallelWWnTransportProcess()
{
  delete fImportancePostStepDoIt;
}

G4VParticleChange *G4ParallelWWnTransportProcess::
PostStepDoIt(const G4Track& aTrack, const G4Step &aStep)
{
  if (aTrack.GetTrackStatus()==fStopAndKill) {
    G4cout << "G4ParallelWWnTransportProcess::PostStepDoIt StopAndKill" << G4endl;
  }
  G4ParallelTransport::PostStepDoIt(aTrack, aStep);

  // get new weight and number of clones
  G4Nsplit_Weight nw(fWeightWindowExaminer.Examine(aTrack.GetWeight(),
						   aTrack.
						   GetKineticEnergy()));

  fImportancePostStepDoIt->DoIt(aTrack, fParticleChange, nw);
  return fParticleChange;
}
  
void G4ParallelWWnTransportProcess::Error(const G4String &m)
{
  G4cout << "ERROR - G4WWnTransportProcess::" << m << G4endl;
  G4Exception("Program aborted.");
}



void G4ParallelWWnTransportProcess::KillTrack() const{
  fParticleChange->SetStatusChange(fStopAndKill);
}


const G4String &G4ParallelWWnTransportProcess::GetName() const {
  return G4ParallelTransport::GetProcessName();
}
