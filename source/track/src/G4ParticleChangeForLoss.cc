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
// $Id: G4ParticleChangeForLoss.cc,v 1.13 2004/06/15 08:17:38 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
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
//
// ------------------------------------------------------------
//
#include "G4ParticleChangeForLoss.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4DynamicParticle.hh"
#include "G4ExceptionSeverity.hh"

G4ParticleChangeForLoss::G4ParticleChangeForLoss():G4VParticleChange()
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

G4ParticleChangeForLoss::G4ParticleChangeForLoss(
             const G4ParticleChangeForLoss &right): G4VParticleChange(right)
{
   if (verboseLevel>1) {
    G4cout << "G4ParticleChangeForLoss::  copy constructor is called " << G4endl;
   }
      currentTrack = right.currentTrack;
      proposedKinEnergy = right.proposedKinEnergy;
      currentCharge = right.currentCharge;
      //theProposedWeight = right.theProposedWeight;
      proposedMomentumDirection = right.proposedMomentumDirection;
}

// assignment operator
G4ParticleChangeForLoss & G4ParticleChangeForLoss::operator=(
                                   const G4ParticleChangeForLoss &right)
{
   if (verboseLevel>1) {
    G4cout << "G4ParticleChangeForLoss:: assignment operator is called " << G4endl;
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
      currentCharge = right.currentCharge;
      //theProposedWeight = right.theProposedWeight;
      proposedMomentumDirection = right.proposedMomentumDirection;
   }
   return *this;
}

//----------------------------------------------------------------
// methods for updating G4Step
//

G4Step* G4ParticleChangeForLoss::UpdateStepForAlongStep(G4Step* pStep)
{
  // A physics process always calculates the final state of the
  // particle relative to the initial state at the beginning
  // of the Step, i.e., based on information of G4Track (or
  // equivalently the PreStepPoint).
  // So, the differences (delta) between these two states have to be
  // calculated and be accumulated in PostStepPoint.

  G4StepPoint* pPreStepPoint = pStep->GetPreStepPoint();
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();

  // calculate new kinetic energy
  G4double kinEnergy = pPostStepPoint->GetKineticEnergy()
                    + (proposedKinEnergy - pPreStepPoint->GetKineticEnergy());

  // update kinetic energy and momentum direction
  if (kinEnergy < 0.0) {
    theLocalEnergyDeposit += kinEnergy;
    kinEnergy = 0.0;
  }

  pPostStepPoint->SetKineticEnergy( kinEnergy );
  pPostStepPoint->SetCharge( currentCharge );

  // update weight 
  // this feature is commented out, it should be overwritten in case
  // if energy loss processes will use biasing
  // G4double newWeight = theProposedWeight/(pPreStepPoint->GetWeight())*(pPostStepPoint->GetWeight());
  // pPostStepPoint->SetWeight( newWeight );


// Not necessary to check now
//#ifdef G4VERBOSE
//  if (debugFlag) CheckIt(*aTrack);
//#endif

  //  Update the G4Step specific attributes
  pStep->AddTotalEnergyDeposit( theLocalEnergyDeposit );
  return pStep;
}

G4Step* G4ParticleChangeForLoss::UpdateStepForPostStep(G4Step* pStep)
{
  // A physics process always calculates the final state of the
  // particle relative to the initial state at the beginning
  // of the Step, i.e., based on information of G4Track (or
  // equivalently the PreStepPoint).
  // So, the differences (delta) between these two states have to be
  // calculated and be accumulated in PostStepPoint.

  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();

  pPostStepPoint->SetKineticEnergy( proposedKinEnergy );
  pPostStepPoint->SetCharge( currentCharge );
  pPostStepPoint->SetMomentumDirection( proposedMomentumDirection );

  // update weight
  // this feature is commented out, it should be overwritten in case
  // if energy loss processes will use biasing
  // pPostStepPoint->SetWeight( theProposedWeight );

// Not necessary to check now
//#ifdef G4VERBOSE
//  if (debugFlag) CheckIt(*aTrack);
//#endif

  //  Update the G4Step specific attributes
  pStep->AddTotalEnergyDeposit( theLocalEnergyDeposit );
  return pStep;
}

//----------------------------------------------------------------
// methods for printing messages
//

void G4ParticleChangeForLoss::DumpInfo() const
{
// use base-class DumpInfo
  G4VParticleChange::DumpInfo();

  G4cout.precision(3);
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
