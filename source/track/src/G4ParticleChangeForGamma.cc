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
// $Id: G4ParticleChangeForGamma.cc 68795 2013-04-05 13:24:46Z gcosmo $
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
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4DynamicParticle.hh"
#include "G4ExceptionSeverity.hh"

G4ParticleChangeForGamma::G4ParticleChangeForGamma()
 : G4VParticleChange(), currentTrack(0), proposedKinEnergy(0.)
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
#ifdef G4VERBOSE
  if (verboseLevel>1) 
    G4cout << "G4ParticleChangeForGamma:: assignment operator is called " << G4endl;
#endif
  
  if (this != &right) {
     if (theNumberOfSecondaries>0) {
#ifdef G4VERBOSE
       if (verboseLevel>0) {
	 G4cout << "G4ParticleChangeForGamma: assignment operator Warning  ";
	 G4cout << "theListOfSecondaries is not empty ";
       }
#endif
       for (G4int index= 0; index<theNumberOfSecondaries; index++){
	 if ( (*theListOfSecondaries)[index] ) delete (*theListOfSecondaries)[index] ;
       }
     }
     delete theListOfSecondaries; 
    theListOfSecondaries =  new G4TrackFastVector();
    theNumberOfSecondaries = right.theNumberOfSecondaries;
    for (G4int index = 0; index<theNumberOfSecondaries; index++){
      G4Track* newTrack =  new G4Track(*((*right.theListOfSecondaries)[index] ));
      theListOfSecondaries->SetElement(index, newTrack);			    }
 
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

void G4ParticleChangeForGamma::AddSecondary(G4DynamicParticle* aParticle)
{
  //  create track
  G4Track* aTrack = new G4Track(aParticle, currentTrack->GetGlobalTime(),
                                           currentTrack->GetPosition());

  //   Touchable handle is copied to keep the pointer
  aTrack->SetTouchableHandle(currentTrack->GetTouchableHandle());

  //  add a secondary
  G4VParticleChange::AddSecondary(aTrack);
}

//----------------------------------------------------------------
// method for updating G4Step
//

G4Step* G4ParticleChangeForGamma::UpdateStepForAtRest(G4Step* pStep)
{
  pStep->AddTotalEnergyDeposit( theLocalEnergyDeposit );
  pStep->SetStepLength( 0.0 );

  if (isParentWeightProposed ){
    pStep->GetPostStepPoint()->SetWeight( theParentWeight );
  }

  return pStep;
}

G4Step* G4ParticleChangeForGamma::UpdateStepForPostStep(G4Step* pStep)
{
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();
  G4Track* pTrack = pStep->GetTrack();

  pPostStepPoint->SetKineticEnergy( proposedKinEnergy );
  pPostStepPoint->SetMomentumDirection( proposedMomentumDirection );
  pPostStepPoint->SetPolarization( proposedPolarization );

  // update velocity for scattering process and particles with mass
  if(proposedKinEnergy > 0.0) { 
    if(pTrack->GetParticleDefinition()->GetPDGMass() > 0.0) {
      pPostStepPoint->SetVelocity(pTrack->CalculateVelocity());
    }
  } 

  if (isParentWeightProposed ){
    pPostStepPoint->SetWeight( theParentWeight );
  }
   
  pStep->AddTotalEnergyDeposit( theLocalEnergyDeposit );
  pStep->AddNonIonizingEnergyDeposit( theNonIonizingEnergyDeposit );
  return pStep;
}

//----------------------------------------------------------------
// methods for printing messages
//

void G4ParticleChangeForGamma::DumpInfo() const
{
// use base-class DumpInfo
  G4VParticleChange::DumpInfo();

  G4int oldprc = G4cout.precision(3);
  G4cout << "        Kinetic Energy (MeV): "
       << std::setw(20) << proposedKinEnergy/MeV
       << G4endl;
  G4cout << "        Momentum Direction: "
       << std::setw(20) << proposedMomentumDirection
       << G4endl;
  G4cout << "        Polarization: "
       << std::setw(20) << proposedPolarization
       << G4endl;
  G4cout.precision(oldprc);
}

G4bool G4ParticleChangeForGamma::CheckIt(const G4Track& aTrack)
{
  G4bool    itsOK = true;
  G4bool    exitWithError = false;

  G4double  accuracy;

  // Energy should not be lager than initial value
  accuracy = ( proposedKinEnergy - aTrack.GetKineticEnergy())/MeV;
  if (accuracy > accuracyForWarning) {
    itsOK = false;
    exitWithError = (accuracy > accuracyForException);
#ifdef G4VERBOSE
    G4cout << "G4ParticleChangeForGamma::CheckIt: ";
    G4cout << "KinEnergy become larger than the initial value!" 
	   << "  Difference:  " << accuracy  << "[MeV] " <<G4endl;
    G4cout << aTrack.GetDefinition()->GetParticleName()
	   << " E=" << aTrack.GetKineticEnergy()/MeV
	   << " pos=" << aTrack.GetPosition().x()/m
	   << ", " << aTrack.GetPosition().y()/m
	   << ", " << aTrack.GetPosition().z()/m
	   << G4endl;
#endif
  }

  // dump out information of this particle change
#ifdef G4VERBOSE
  if (!itsOK) DumpInfo();
#endif

  // Exit with error
  if (exitWithError) {
    G4Exception("G4ParticleChangeForGamma::CheckIt",
		"TRACK004",EventMustBeAborted,
		"energy was  illegal");
  }

  //correction
  if (!itsOK) {
    proposedKinEnergy = aTrack.GetKineticEnergy();
  }

  itsOK = (itsOK) && G4VParticleChange::CheckIt(aTrack);
  return itsOK;
}
