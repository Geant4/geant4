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
// $Id: G4ParticleChangeForTransport.cc 87698 2014-12-17 09:41:28Z gcosmo $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	
//
// ------------------------------------------------------------
//   Implemented for the new scheme                 10 May. 1998  H.Kurahige
//   Correct tratment of fpNextTouchable            12 May. 1998  H.Kurashige
//   Change to the next volume only if energy>0     19 Jan. 2004  V.Ivanchenko
// --------------------------------------------------------------

#include "G4ParticleChangeForTransport.hh"
#include "G4TouchableHandle.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4TrackFastVector.hh"
#include "G4DynamicParticle.hh"

G4ParticleChangeForTransport::G4ParticleChangeForTransport()
  : G4ParticleChange(), isMomentumChanged(false), theMaterialChange(0),
    theMaterialCutsCoupleChange(0), theSensitiveDetectorChange(0),
    fpVectorOfAuxiliaryPointsPointer(0)
{
  if (verboseLevel>2) {
    G4cout << "G4ParticleChangeForTransport::G4ParticleChangeForTransport() "
           << G4endl;
  }
}

G4ParticleChangeForTransport::~G4ParticleChangeForTransport() 
{
  if (verboseLevel>2) {
    G4cout << "G4ParticleChangeForTransport::~G4ParticleChangeForTransport() "
           << G4endl;
  }
}

G4ParticleChangeForTransport::
G4ParticleChangeForTransport(const G4ParticleChangeForTransport &r)
  : G4ParticleChange(r),
    fpVectorOfAuxiliaryPointsPointer(0)
{
  if (verboseLevel>0) {
    G4cout << "G4ParticleChangeForTransport::  copy constructor is called "
           << G4endl;
  }
  theTouchableHandle = r.theTouchableHandle;
  isMomentumChanged = r.isMomentumChanged;
  theMaterialChange = r.theMaterialChange;
  theMaterialCutsCoupleChange = r.theMaterialCutsCoupleChange;
  theSensitiveDetectorChange = r.theSensitiveDetectorChange;
}

// assignemnt operator
G4ParticleChangeForTransport &
G4ParticleChangeForTransport::operator=(const G4ParticleChangeForTransport &r)
{
   if (verboseLevel>1) {
    G4cout << "G4ParticleChangeForTransport:: assignment operator is called "
           << G4endl;
   }
   if (this != &r)
   {
      theListOfSecondaries = r.theListOfSecondaries;
      theSizeOftheListOfSecondaries = r.theSizeOftheListOfSecondaries;
      theNumberOfSecondaries = r.theNumberOfSecondaries;
      theStatusChange = r.theStatusChange;
      theTouchableHandle = r.theTouchableHandle;
      theMaterialChange = r.theMaterialChange;
      theMaterialCutsCoupleChange = r.theMaterialCutsCoupleChange;
      theSensitiveDetectorChange = r.theSensitiveDetectorChange;
      theMomentumDirectionChange = r.theMomentumDirectionChange;
      thePolarizationChange = r.thePolarizationChange;
      thePositionChange = r.thePositionChange;
      theTimeChange = r.theTimeChange;
      theEnergyChange = r.theEnergyChange;
      theVelocityChange        = r.theVelocityChange;
      theTrueStepLength = r.theTrueStepLength;
      theLocalEnergyDeposit = r.theLocalEnergyDeposit;
      theSteppingControlFlag = r.theSteppingControlFlag;
   }
   return *this;
}

//----------------------------------------------------------------
// methods for updating G4Step
//

G4Step* G4ParticleChangeForTransport::UpdateStepForAtRest(G4Step* pStep)
{
  // Nothing happens for AtRestDoIt
  if (verboseLevel>0) {
    G4cout << "G4ParticleChangeForTransport::UpdateStepForAtRest() is called"
           << G4endl;
    G4cout << " Nothing happens for this method " << G4endl;
  }
  //  Update the G4Step specific attributes
  return UpdateStepInfo(pStep);
}


G4Step* G4ParticleChangeForTransport::UpdateStepForAlongStep(G4Step* pStep)
{
  // Smooth curved tajectory representation: let the Step know about
  // the auxiliary trajectory points (jacek 30/10/2002)
  pStep->SetPointerToVectorOfAuxiliaryPoints(fpVectorOfAuxiliaryPointsPointer);

  // copy of G4ParticleChange::UpdateStepForAlongStep
  //  i.e. no effect for touchable

  // A physics process always calculates the final state of the
  // particle relative to the initial state at the beginning
  // of the Step, i.e., based on information of G4Track (or
  // equivalently the PreStepPoint).
  // So, the differences (delta) between these two states have to be
  // calculated and be accumulated in PostStepPoint.

  // Take note that the return type of GetMomentumChange is a
  // pointer to G4ThreeVector. Also it is a normalized
  // momentum vector.

  G4StepPoint* pPreStepPoint  = pStep->GetPreStepPoint();
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();
  G4Track*     aTrack  = pStep->GetTrack();
  G4double     mass = aTrack->GetDynamicParticle()->GetMass();

  // uodate kinetic energy
  //  now assume that no energy change in transportation
  //  However it is not true in electric fields
  //  Case for changing energy will be implemented in future


  // update momentum direction and energy
  if (isMomentumChanged) {
    G4double energy;
    energy= pPostStepPoint->GetKineticEnergy()
                 + (theEnergyChange - pPreStepPoint->GetKineticEnergy());

    // calculate new momentum
    G4ThreeVector pMomentum =  pPostStepPoint->GetMomentum()
                     + ( CalcMomentum(theEnergyChange, theMomentumDirectionChange, mass)
	                  - pPreStepPoint->GetMomentum());
    G4double      tMomentum = pMomentum.mag();
    G4ThreeVector direction(1.0,0.0,0.0);
    if( tMomentum > 0. ){
      G4double  inv_Momentum= 1.0 / tMomentum;
      direction= pMomentum * inv_Momentum;
    }
    pPostStepPoint->SetMomentumDirection(direction);
    pPostStepPoint->SetKineticEnergy( energy );
  }
  if (isVelocityChanged)  pPostStepPoint->SetVelocity(theVelocityChange);

  // stop case should not occur
  //pPostStepPoint->SetMomentumDirection(G4ThreeVector(1., 0., 0.));


  // update polarization
  pPostStepPoint->AddPolarization( thePolarizationChange
  				   - pPreStepPoint->GetPolarization());

  // update position and time
  pPostStepPoint->AddPosition( thePositionChange
			       - pPreStepPoint->GetPosition() );
  pPostStepPoint->AddGlobalTime( theTimeChange
				 - pPreStepPoint->GetLocalTime());
  pPostStepPoint->AddLocalTime( theTimeChange
				 - pPreStepPoint->GetLocalTime());
  pPostStepPoint->AddProperTime( theProperTimeChange
				 - pPreStepPoint->GetProperTime());

#ifdef G4VERBOSE
  if (debugFlag) CheckIt(*aTrack);
#endif

  //  Update the G4Step specific attributes
  //pStep->SetStepLength( theTrueStepLength );
  //  pStep->AddTotalEnergyDeposit( theLocalEnergyDeposit );
  pStep->SetControlFlag( theSteppingControlFlag );
  return pStep;
  //  return UpdateStepInfo(pStep);
}

G4Step* G4ParticleChangeForTransport::UpdateStepForPostStep(G4Step* pStep)
{
  // A physics process always calculates the final state of the particle

  // Change volume only if some kinetic energy remains
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();
  if(pPostStepPoint->GetKineticEnergy() > 0.0) {

    // update next touchable
    // (touchable can be changed only at PostStepDoIt)
    pPostStepPoint->SetTouchableHandle( theTouchableHandle );

    pPostStepPoint->SetMaterial( theMaterialChange );
    pPostStepPoint->SetMaterialCutsCouple( theMaterialCutsCoupleChange );
    pPostStepPoint->SetSensitiveDetector( theSensitiveDetectorChange );
  }
  if( this->GetFirstStepInVolume() ){
    pStep->SetFirstStepFlag();
  }else{
    pStep->ClearFirstStepFlag(); 
  }  
  if( this->GetLastStepInVolume() ){
    pStep->SetLastStepFlag();
  }else{
    pStep->ClearLastStepFlag(); 
  }
  // It used to call base class's method
  //   - but this would copy uninitialised data members
  // return G4ParticleChange::UpdateStepForPostStep(pStep);

  // Copying what the base class does would instead
  //   - also not useful
  // return G4VParticleChange::UpdateStepInfo(pStep);

  return pStep;
}


//----------------------------------------------------------------
// methods for printing messages
//

void G4ParticleChangeForTransport::DumpInfo() const
{
// use base-class DumpInfo
  G4ParticleChange::DumpInfo();

  G4int oldprc = G4cout.precision(3);
  G4cout << "        Touchable (pointer) : " 
         << std::setw(20) << theTouchableHandle() << G4endl; 
  G4cout.precision(oldprc);
}

