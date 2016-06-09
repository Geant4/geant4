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
// $Id: G4ParallelWWnTransportProcess.cc,v 1.6 2006/06/29 21:12:18 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ParallelWWnTransportProcess.cc
//
// ----------------------------------------------------------------------

#include "G4Types.hh"
#include "G4ParallelWWnTransportProcess.hh"
#include "G4VWeightWindowExaminer.hh"
#include "G4VTrackTerminator.hh"
#include "G4SamplingPostStepAction.hh"

G4ParallelWWnTransportProcess::G4ParallelWWnTransportProcess(
                      const G4VWeightWindowExaminer &aWeightWindowExaminer,
                            G4VPGeoDriver &pgeodriver, 
                            G4VParallelStepper &aStepper,
                      const G4VTrackTerminator *TrackTerminator,
                      const G4String &aName ) 
  : G4ParallelTransport(pgeodriver, aStepper, aName),
    fParticleChange(G4ParallelTransport::fParticleChange),
    fWeightWindowExaminer(aWeightWindowExaminer),
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

G4ParallelWWnTransportProcess::~G4ParallelWWnTransportProcess()
{
  delete fPostStepAction;
}

G4VParticleChange *
G4ParallelWWnTransportProcess::PostStepDoIt(const G4Track& aTrack,
                                            const G4Step &aStep)
{
  if (aTrack.GetTrackStatus()==fStopAndKill)
  {
    G4cout << "WARNING - G4ParallelWWnTransportProcess::PostStepDoIt()"
           << "          StopAndKill track !"
           << G4endl;
  }
  G4ParallelTransport::PostStepDoIt(aTrack, aStep);

  // get new weight and number of clones
  G4Nsplit_Weight nw(fWeightWindowExaminer.Examine(aTrack.GetWeight(),
                                                   aTrack.
                                                   GetKineticEnergy()));

  fPostStepAction->DoIt(aTrack, fParticleChange, nw);
  return fParticleChange;
}
  
void G4ParallelWWnTransportProcess::Error(const G4String &m)
{
  G4Exception("G4ParallelWWnTransportProcess::Error()",
              "ProgramError", FatalException, m);
}

void G4ParallelWWnTransportProcess::KillTrack() const
{
  fParticleChange->ProposeTrackStatus(fStopAndKill);
}

const G4String &G4ParallelWWnTransportProcess::GetName() const
{
  return G4ParallelTransport::GetProcessName();
}
