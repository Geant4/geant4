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
// $Id: G4ParticleChangeForDecay.cc,v 1.12 2010-07-21 09:30:15 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	
//	
// ------------------------------------------------------------
//   Implemented for the new scheme                 23 Mar. 1998  H.Kurahige
//   Remove modification of energy/momentum         20 Jul, 1998  H.Kurashige
// --------------------------------------------------------------

#include "G4ParticleChangeForDecay.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4TrackFastVector.hh"
#include "G4DynamicParticle.hh"
#include "G4ExceptionSeverity.hh"

G4ParticleChangeForDecay::G4ParticleChangeForDecay()
  : G4VParticleChange(), theTimeChange(0.)
{
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4ParticleChangeForDecay::G4ParticleChangeForDecay() " << G4endl;
  }
#endif
}

G4ParticleChangeForDecay::~G4ParticleChangeForDecay() 
{
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4ParticleChangeForDecay::~G4ParticleChangeForDecay() " << G4endl;
  }
#endif
}

// copy and assignment operators are implemented as "shallow copy"
G4ParticleChangeForDecay::G4ParticleChangeForDecay(const G4ParticleChangeForDecay &right): G4VParticleChange(right)
{
  theTimeChange = right.theTimeChange;
  thePolarizationChange = right.thePolarizationChange;
}


G4ParticleChangeForDecay & G4ParticleChangeForDecay::operator=(const G4ParticleChangeForDecay &right)
{
   if (this != &right){
      theListOfSecondaries = right.theListOfSecondaries;
      theSizeOftheListOfSecondaries = right.theSizeOftheListOfSecondaries;
      theNumberOfSecondaries = right.theNumberOfSecondaries;
      theStatusChange = right.theStatusChange;
      theTrueStepLength = right.theTrueStepLength;
      theLocalEnergyDeposit = right.theLocalEnergyDeposit;
      theSteppingControlFlag = right.theSteppingControlFlag;
     
      theTimeChange = right.theTimeChange;
      thePolarizationChange = right.thePolarizationChange;
   }
   return *this;
}

G4bool G4ParticleChangeForDecay::operator==(const G4ParticleChangeForDecay &right) const
{
   return ((G4VParticleChange *)this == (G4VParticleChange *) &right);
}

G4bool G4ParticleChangeForDecay::operator!=(const G4ParticleChangeForDecay &right) const
{
   return ((G4VParticleChange *)this != (G4VParticleChange *) &right);
}

//----------------------------------------------------------------
// methods for Initialization
//
void G4ParticleChangeForDecay::Initialize(const G4Track& track)
{
  // use base class's method at first
  G4VParticleChange::Initialize(track);

  const G4DynamicParticle*  pParticle = track.GetDynamicParticle();

  // set Time e equal to those of the parent track
  theTimeChange          = track.GetGlobalTime();

  // set the Polarization equal to those of the parent track
  thePolarizationChange  = pParticle->GetPolarization();
}

//----------------------------------------------------------------
// methods for updating G4Step 
//

G4Step* G4ParticleChangeForDecay::UpdateStepForPostStep(G4Step* pStep)
{ 
  //  Update the G4Step specific attributes 
  return UpdateStepInfo(pStep);
}


G4Step* G4ParticleChangeForDecay::UpdateStepForAtRest(G4Step* pStep)
{ 
  // A physics process always calculates the final state of the particle

  G4StepPoint* pPreStepPoint  = pStep->GetPreStepPoint(); 
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint(); 

  // update polarization
  pPostStepPoint->SetPolarization( thePolarizationChange );
 
  // update time
  pPostStepPoint->SetGlobalTime( theTimeChange  );
  pPostStepPoint->AddLocalTime(  theTimeChange 
                                 - pPreStepPoint->GetGlobalTime());
  pPostStepPoint->AddProperTime( theTimeChange  
                                 - pPreStepPoint->GetGlobalTime());


#ifdef G4VERBOSE
  G4Track*     aTrack  = pStep->GetTrack();
  if (debugFlag) CheckIt(*aTrack);
#endif

  //  Update the G4Step specific attributes 
  return UpdateStepInfo(pStep);
}

void G4ParticleChangeForDecay::DumpInfo() const
{
// Show header
  G4VParticleChange::DumpInfo();

  G4int oldprc = G4cout.precision(3);
  G4cout << "        Time (ns)           : " 
         << std::setw(20) << theTimeChange/ns << G4endl;
  G4cout.precision(oldprc);
}

G4bool G4ParticleChangeForDecay::CheckIt(const G4Track& aTrack)
{
  G4bool    exitWithError = false;

  G4double  accuracy;

  // global time should not go back
  G4bool itsOK =true;
  accuracy = -1.0*(theTimeChange - aTrack.GetGlobalTime())/ns;
  if (accuracy > accuracyForWarning) {
#ifdef G4VERBOSE
    G4cout << "  G4ParticleChangeForDecay::CheckIt    : ";
    G4cout << "the global time goes back  !!" << G4endl;
    G4cout << "  Difference:  " << accuracy  << "[ns] " <<G4endl;
#endif
    itsOK = false;
    if (accuracy > accuracyForException) exitWithError = true;
  }

  // dump out information of this particle change
#ifdef G4VERBOSE
  if (!itsOK) { 
    G4cout << " G4ParticleChangeForDecay::CheckIt " <<G4endl;
    DumpInfo();
  }
#endif

  // Exit with error
  if (exitWithError) {
    G4Exception("G4ParticleChangeForDecay::CheckIt",
		"500",
		EventMustBeAborted,
		"time was  illegal");
  } 

  // correction
  if (!itsOK) {
    theTimeChange = aTrack.GetGlobalTime();
  }

  itsOK = (itsOK) && G4VParticleChange::CheckIt(aTrack);
  return itsOK;
}





