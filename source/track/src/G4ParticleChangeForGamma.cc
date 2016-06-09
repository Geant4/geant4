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
// $Id: G4ParticleChangeForGamma.cc,v 1.3 2006/08/28 16:10:29 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
//
// --------------------------------------------------------------
//	GEANT 4 class implementation file
//
// ------------------------------------------------------------
//   15 April 2005 V.Ivanchenko for gamma EM processes
// --------------------------------------------------------------
//
//   Modified::
//   28.08.06 V.Ivanchenko Add access to current track and polarizaion
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
  if (verboseLevel>1) 
    G4cout << "G4ParticleChangeForGamma::  copy constructor is called " << G4endl;
  
  currentTrack = right.currentTrack;
  proposedKinEnergy = right.proposedKinEnergy;
  proposedMomentumDirection = right.proposedMomentumDirection;
  proposedPolarization = right.proposedPolarization;
}

// assignment operator
G4ParticleChangeForGamma & G4ParticleChangeForGamma::operator=(
			   const G4ParticleChangeForGamma &right)
{
  if (verboseLevel>1) 
    G4cout << "G4ParticleChangeForGamma:: assignment operator is called " << G4endl;
  
  if (this != &right) {
    theListOfSecondaries = right.theListOfSecondaries;
    theSizeOftheListOfSecondaries = right.theSizeOftheListOfSecondaries;
    theNumberOfSecondaries = right.theNumberOfSecondaries;
    theStatusChange = right.theStatusChange;
    theLocalEnergyDeposit = right.theLocalEnergyDeposit;
    theSteppingControlFlag = right.theSteppingControlFlag;
    theParentWeight = right.theParentWeight;

    currentTrack = right.currentTrack;
    proposedKinEnergy = right.proposedKinEnergy;
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
