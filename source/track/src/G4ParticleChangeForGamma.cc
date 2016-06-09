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
// $Id: G4ParticleChangeForGamma.cc,v 1.1 2005/04/14 18:00:15 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//
// --------------------------------------------------------------
//	GEANT 4 class implementation file
//
// ------------------------------------------------------------
//   15 April 2005 V.Ivanchenko for gamma EM processes
// --------------------------------------------------------------
//
//   Modified:
//
// ------------------------------------------------------------
//
#include "G4ParticleChangeForGamma.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4DynamicParticle.hh"
#include "G4ExceptionSeverity.hh"

G4ParticleChangeForGamma::G4ParticleChangeForGamma():G4VParticleChange()
{
  theSteppingControlFlag = NormalCondition;
  debugFlag = false;
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4ParticleChangeForGamma::G4ParticleChangeForGamma() " << G4endl;
  }
#endif
}

G4ParticleChangeForGamma::~G4ParticleChangeForGamma()
{
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4ParticleChangeForGamma::~G4ParticleChangeForGamma() " << G4endl;
  }
#endif
}

G4ParticleChangeForGamma::G4ParticleChangeForGamma(
             const G4ParticleChangeForGamma &right): G4VParticleChange(right)
{
   if (verboseLevel>1) {
    G4cout << "G4ParticleChangeForGamma::  copy constructor is called " << G4endl;
   }
      currentTrack = right.currentTrack;
      proposedKinEnergy = right.proposedKinEnergy;
      //theProposedWeight = right.theProposedWeight;
      proposedMomentumDirection = right.proposedMomentumDirection;
      proposedPolarization = right.proposedPolarization;
}

// assignment operator
G4ParticleChangeForGamma & G4ParticleChangeForGamma::operator=(
                                   const G4ParticleChangeForGamma &right)
{
   if (verboseLevel>1) {
    G4cout << "G4ParticleChangeForGamma:: assignment operator is called " << G4endl;
   }
   if (this != &right)
   {
      theListOfSecondaries = right.theListOfSecondaries;
      theSizeOftheListOfSecondaries = right.theSizeOftheListOfSecondaries;
      theNumberOfSecondaries = right.theNumberOfSecondaries;
      theStatusChange = right.theStatusChange;
      theLocalEnergyDeposit = right.theLocalEnergyDeposit;
      theSteppingControlFlag = right.theSteppingControlFlag;

      currentTrack = right.currentTrack;
      proposedKinEnergy = right.proposedKinEnergy;
      //theProposedWeight = right.theProposedWeight;
      proposedMomentumDirection = right.proposedMomentumDirection;
      proposedPolarization = right.proposedPolarization;
   }
   return *this;
}

//----------------------------------------------------------------
// methods for printing messages
//

void G4ParticleChangeForGamma::DumpInfo() const
{
// use base-class DumpInfo
  G4VParticleChange::DumpInfo();

  G4cout.precision(3);
  G4cout << "        Kinetic Energy (MeV): "
       << std::setw(20) << proposedKinEnergy/MeV
       << G4endl;
  G4cout << "        Momentum Direction: "
       << std::setw(20) << proposedMomentumDirection
       << G4endl;
  G4cout << "        Polarization: "
       << std::setw(20) << proposedPolarization
       << G4endl;
}

G4bool G4ParticleChangeForGamma::CheckIt(const G4Track& aTrack)
{
  G4bool    itsOK = true;
  G4bool    exitWithError = false;

  G4double  accuracy;

  // Energy should not be lager than initial value
  accuracy = ( proposedKinEnergy - aTrack.GetKineticEnergy())/MeV;
  if (accuracy > accuracyForWarning) {
#ifdef G4VERBOSE
    G4cout << "G4ParticleChangeForGamma::CheckIt: ";
    G4cout << "KinEnergy become larger than the initial value!" << G4endl;
    G4cout << "  Difference:  " << accuracy  << "[MeV] " <<G4endl;
#endif
    itsOK = false;
    if (accuracy > accuracyForException) exitWithError = true;
  }

  // dump out information of this particle change
#ifdef G4VERBOSE
  if (!itsOK) {
    G4cout << "G4ParticleChangeForGamma::CheckIt " << G4endl;
    DumpInfo();
  }
#endif

  // Exit with error
  if (exitWithError) {
    G4Exception("G4ParticleChangeForGamma::CheckIt",
		"400",
		EventMustBeAborted,
		"energy was  illegal");
  }

  //correction
  if (!itsOK) {
    proposedKinEnergy = aTrack.GetKineticEnergy();
  }

  itsOK = (itsOK) && G4VParticleChange::CheckIt(aTrack);
  return itsOK;
}
