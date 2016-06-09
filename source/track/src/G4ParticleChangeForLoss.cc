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
// $Id: G4ParticleChangeForLoss.cc,v 1.18 2010-07-21 09:30:15 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// --------------------------------------------------------------
//	GEANT 4 class implementation file
//
// ------------------------------------------------------------
//   Implemented for the new scheme                 23 Mar. 1998  H.Kurahige
// --------------------------------------------------------------
//
//   Modified:
//   16.01.04 V.Ivanchenko update for model variant of energy loss
//   15.04.05 V.Ivanchenko inline update methods
//   28.08.06 V.Ivanchenko Add access to current track and polarizaion
//
// ------------------------------------------------------------
//
#include "G4ParticleChangeForLoss.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4DynamicParticle.hh"
#include "G4ExceptionSeverity.hh"

G4ParticleChangeForLoss::G4ParticleChangeForLoss()
  : G4VParticleChange(), currentTrack(0), proposedKinEnergy(0.),
    lowEnergyLimit(1.0*eV), currentCharge(0.)
{
  theSteppingControlFlag = NormalCondition;
  debugFlag = false;
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4ParticleChangeForLoss::G4ParticleChangeForLoss() " << G4endl;
  }
#endif
}

G4ParticleChangeForLoss::~G4ParticleChangeForLoss()
{
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4ParticleChangeForLoss::~G4ParticleChangeForLoss() " << G4endl;
  }
#endif
}

G4ParticleChangeForLoss::
G4ParticleChangeForLoss(const G4ParticleChangeForLoss &right)
  : G4VParticleChange(right)
{
  if (verboseLevel>1) {
    G4cout << "G4ParticleChangeForLoss::  copy constructor is called " << G4endl;
  }
  currentTrack = right.currentTrack;
  proposedKinEnergy = right.proposedKinEnergy;
  lowEnergyLimit = right.lowEnergyLimit;
  currentCharge = right.currentCharge;
  proposedMomentumDirection = right.proposedMomentumDirection;
}

// assignment operator
G4ParticleChangeForLoss & G4ParticleChangeForLoss::operator=(
                                   const G4ParticleChangeForLoss &right)
{
  if (verboseLevel>1) {
    G4cout << "G4ParticleChangeForLoss:: assignment operator is called " << G4endl;
  }
  if (this != &right) {
    theListOfSecondaries = right.theListOfSecondaries;
    theSizeOftheListOfSecondaries = right.theSizeOftheListOfSecondaries;
    theNumberOfSecondaries = right.theNumberOfSecondaries;
    theStatusChange = right.theStatusChange;
    theLocalEnergyDeposit = right.theLocalEnergyDeposit;
    theSteppingControlFlag = right.theSteppingControlFlag;
    theParentWeight = right.theParentWeight;
    isParentWeightProposed = right.isParentWeightProposed;
    isParentWeightSetByProcess = right.isParentWeightSetByProcess;
    fSetSecondaryWeightByProcess = right.fSetSecondaryWeightByProcess;

    currentTrack = right.currentTrack;
    proposedKinEnergy = right.proposedKinEnergy;
    currentCharge = right.currentCharge;
    proposedMomentumDirection = right.proposedMomentumDirection;
  }
  return *this;
}

//----------------------------------------------------------------
// methods for printing messages
//

void G4ParticleChangeForLoss::DumpInfo() const
{
// use base-class DumpInfo
  G4VParticleChange::DumpInfo();

  G4int oldprc = G4cout.precision(3);
  G4cout << "        Charge (eplus)   : "
       << std::setw(20) << currentCharge/eplus
       << G4endl;
  G4cout << "        Kinetic Energy (MeV): "
       << std::setw(20) << proposedKinEnergy/MeV
       << G4endl;
  G4cout << "        Momentum Direct - x : "
       << std::setw(20) << proposedMomentumDirection.x()
       << G4endl;
  G4cout << "        Momentum Direct - y : "
       << std::setw(20) << proposedMomentumDirection.y()
       << G4endl;
  G4cout << "        Momentum Direct - z : "
       << std::setw(20) << proposedMomentumDirection.z()
       << G4endl;
  G4cout.precision(oldprc);
}

G4bool G4ParticleChangeForLoss::CheckIt(const G4Track& aTrack)
{
  G4bool    itsOK = true;
  G4bool    exitWithError = false;

  G4double  accuracy;

  // Energy should not be lager than initial value
  accuracy = ( proposedKinEnergy - aTrack.GetKineticEnergy())/MeV;
  if (accuracy > accuracyForWarning) {
#ifdef G4VERBOSE
    G4cout << "G4ParticleChangeForLoss::CheckIt: ";
    G4cout << "KinEnergy become larger than the initial value!" << G4endl;
    G4cout << "  Difference:  " << accuracy  << "[MeV] " <<G4endl;
#endif
    itsOK = false;
    if (accuracy > accuracyForException) exitWithError = true;
  }

  // dump out information of this particle change
#ifdef G4VERBOSE
  if (!itsOK) {
    G4cout << "G4ParticleChangeForLoss::CheckIt " << G4endl;
    DumpInfo();
  }
#endif

  // Exit with error
  if (exitWithError) {
    G4Exception("G4ParticleChangeForLoss::CheckIt",
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

//----------------------------------------------------------------
// methods for updating G4Step
//

G4Step* G4ParticleChangeForLoss::UpdateStepForAlongStep(G4Step* pStep)
{
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();

  // accumulate change of the kinetic energy
  G4double kinEnergy = pPostStepPoint->GetKineticEnergy() +
    (proposedKinEnergy - pStep->GetPreStepPoint()->GetKineticEnergy());

  // update kinetic energy and charge
  if (kinEnergy < lowEnergyLimit) {
    theLocalEnergyDeposit += kinEnergy;
    kinEnergy = 0.0;
  } else {
    pPostStepPoint->SetCharge( currentCharge );
  }
  pPostStepPoint->SetKineticEnergy( kinEnergy );

  if (isParentWeightProposed) {
    // update weight 
    G4StepPoint* pPreStepPoint = pStep->GetPreStepPoint();
    G4double newWeight= theParentWeight/(pPreStepPoint->GetWeight())
      * (pPostStepPoint->GetWeight());
    if (isParentWeightSetByProcess) pPostStepPoint->SetWeight( newWeight );
    if (!fSetSecondaryWeightByProcess) {    
      //Set weight of secondary tracks
      for (G4int index= 0; index<theNumberOfSecondaries; index++){
        if ( (*theListOfSecondaries)[index] ) {
          ((*theListOfSecondaries)[index])->SetWeight(newWeight); ;
        }
      }
    }
  }

  pStep->AddTotalEnergyDeposit( theLocalEnergyDeposit );
  pStep->AddNonIonizingEnergyDeposit( theNonIonizingEnergyDeposit );
  return pStep;
}

G4Step* G4ParticleChangeForLoss::UpdateStepForPostStep(G4Step* pStep)
{
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();
  pPostStepPoint->SetCharge( currentCharge );
  pPostStepPoint->SetMomentumDirection( proposedMomentumDirection );
  pPostStepPoint->SetKineticEnergy( proposedKinEnergy );
  pPostStepPoint->SetPolarization( proposedPolarization );

  if (isParentWeightProposed) {
    if (isParentWeightSetByProcess) pPostStepPoint->SetWeight( theParentWeight );
    if (!fSetSecondaryWeightByProcess) {    
      // Set weight of secondary tracks
      for (G4int index= 0; index<theNumberOfSecondaries; index++){
        if ( (*theListOfSecondaries)[index] ) {
          ((*theListOfSecondaries)[index])->SetWeight( theParentWeight );
        }
      }
    }
  }

  pStep->AddTotalEnergyDeposit( theLocalEnergyDeposit );
  pStep->AddNonIonizingEnergyDeposit( theNonIonizingEnergyDeposit );
  return pStep;
}

void G4ParticleChangeForLoss::AddSecondary(G4DynamicParticle* aParticle)
{
  //  create track
  G4Track* aTrack = new G4Track(aParticle, currentTrack->GetGlobalTime(),
                                           currentTrack->GetPosition());

  //   Touchable handle is copied to keep the pointer
  aTrack->SetTouchableHandle(currentTrack->GetTouchableHandle());

  //  add a secondary
  G4VParticleChange::AddSecondary(aTrack);
}
