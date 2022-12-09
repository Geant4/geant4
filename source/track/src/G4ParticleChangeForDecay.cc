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
// G4ParticleChangeForDecay class implementation
//
// Author: Hisaya Kurashige, 23 March 1998  
// --------------------------------------------------------------------

#include "G4ParticleChangeForDecay.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4DynamicParticle.hh"
#include "G4ExceptionSeverity.hh"

// --------------------------------------------------------------------
G4ParticleChangeForDecay::G4ParticleChangeForDecay()
{}

// --------------------------------------------------------------------
void G4ParticleChangeForDecay::Initialize(const G4Track& track)
{
  // use base class's method at first
  G4VParticleChange::Initialize(track);

  // set TimeChange equal to local time of the parent track
  theLocalTime0 = theTimeChange = track.GetLocalTime();

  // set initial Local/Global time of the parent track
  theGlobalTime0 = track.GetGlobalTime();

  // set the Polarization equal to those of the parent track
  thePolarizationChange = track.GetDynamicParticle()->GetPolarization();
}

// --------------------------------------------------------------------
G4Step* G4ParticleChangeForDecay::UpdateStepForPostStep(G4Step* pStep)
{
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();

  if(isParentWeightProposed)
  {
    pPostStepPoint->SetWeight(theParentWeight);
  }
  // update polarization
  pPostStepPoint->SetPolarization(thePolarizationChange);

  // update the G4Step specific attributes
  return UpdateStepInfo(pStep);
}

// --------------------------------------------------------------------
G4Step* G4ParticleChangeForDecay::UpdateStepForAtRest(G4Step* pStep)
{
  // A physics process always calculates the final state of the particle

  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();

  // update polarization
  pPostStepPoint->SetPolarization(thePolarizationChange);

  // update time
  pPostStepPoint->SetGlobalTime(GetGlobalTime());
  pPostStepPoint->SetLocalTime(theTimeChange);
  pPostStepPoint->AddProperTime(theTimeChange - theLocalTime0);

#ifdef G4VERBOSE
  if(debugFlag) { CheckIt(*theCurrentTrack); }
#endif

  if(isParentWeightProposed)
    pPostStepPoint->SetWeight(theParentWeight);

  //  Update the G4Step specific attributes
  return UpdateStepInfo(pStep);
}

// --------------------------------------------------------------------
void G4ParticleChangeForDecay::DumpInfo() const
{
  // Show header
  G4VParticleChange::DumpInfo();

  G4long oldprc = G4cout.precision(8);
  G4cout << "      -----------------------------------------------" << G4endl;
  G4cout << "    G4ParticleChangeForDecay proposes: " << G4endl;
  G4cout << "    Proposed local Time (ns): " << std::setw(20)
         << theTimeChange / ns << G4endl;
  G4cout << "    Initial local Time (ns) : " << std::setw(20)
         << theLocalTime0 / ns << G4endl;
  G4cout << "    Initial global Time (ns): " << std::setw(20)
         << theGlobalTime0 / ns << G4endl;
  G4cout << "    Current global Time (ns): " << std::setw(20)
         << theCurrentTrack->GetGlobalTime() / ns << G4endl;
  G4cout.precision(oldprc);
}

// --------------------------------------------------------------------
G4bool G4ParticleChangeForDecay::CheckIt(const G4Track& aTrack)
{
  // local time should not go back
  G4bool itsOK = true;
  if(theTimeChange < theLocalTime0)
  {
    itsOK         = false;
    ++nError;
#ifdef G4VERBOSE
    if(nError < maxError)
    {
      G4cout << "  G4ParticleChangeForDecay::CheckIt    : ";
      G4cout << "the local time goes back  !!"
	     << "  Difference:  " << (theTimeChange - theLocalTime0) / ns 
             << "[ns] " << G4endl;
      G4cout << "initial local time " << theLocalTime0 / ns << "[ns] "
	     << "initial global time " << theGlobalTime0 / ns << "[ns] "
	     << G4endl;
    }
#endif
    theTimeChange = theLocalTime0;
  }

  if(!itsOK)
  {
    if(nError < maxError)
    {
#ifdef G4VERBOSE
      // dump out information of this particle change
      DumpInfo();
#endif
      G4Exception("G4ParticleChangeForDecay::CheckIt()", "TRACK005",
                  JustWarning, "time is illegal");
    }
  }

  return (itsOK) && G4VParticleChange::CheckIt(aTrack);
}
