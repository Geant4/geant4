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
// $Id: G4ParticleChange.cc,v 1.32 2010-07-21 09:30:15 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	
//	
// ------------------------------------------------------------
//   Implemented for the new scheme                 23 Mar. 1998  H.Kurahige
//   Change default debug flag to false             10 May. 1998  H.Kurahige
//   Add Track weight                               12 Nov. 1998  H.Kurashige
//   Activate CheckIt method for VERBOSE mode       14 Dec. 1998 H.Kurashige
//   Modified CheckIt method for time                9 Feb. 1999 H.Kurashige
//   Rename SetXXX methods to ProposeXXX   DynamicCharge  Oct. 2005 H.Kurashige
//   Add get/ProposeMagneticMoment                  Mar 2007 H.Kurashige
// --------------------------------------------------------------

#include "G4ParticleChange.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4TrackFastVector.hh"
#include "G4DynamicParticle.hh"
#include "G4ExceptionSeverity.hh"

G4ParticleChange::G4ParticleChange()
  : G4VParticleChange(), theEnergyChange(0.), theTimeChange(0.),
    theProperTimeChange(0.), theMassChange(0.), theChargeChange(0.),
    theMagneticMomentChange(0.), theCurrentTrack(0)
{
}

G4ParticleChange::~G4ParticleChange() 
{
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4ParticleChange::~G4ParticleChange() " << G4endl;
  }
#endif
}

// copy constructor
G4ParticleChange::G4ParticleChange(const G4ParticleChange &right): G4VParticleChange(right)
{
   if (verboseLevel>1) {
    G4cout << "G4ParticleChange::  copy constructor is called " << G4endl;
   }
   theMomentumDirectionChange = right.theMomentumDirectionChange;
   thePolarizationChange = right.thePolarizationChange;
   thePositionChange = right.thePositionChange;
   theTimeChange = right.theTimeChange;
   theEnergyChange = right.theEnergyChange;
   theMassChange = right.theMassChange;
   theChargeChange = right.theChargeChange;
   theMagneticMomentChange = right.theMagneticMomentChange;
   //   theWeightChange = right.theWeightChange;
   theProperTimeChange = right.theProperTimeChange;
}

// assignemnt operator
G4ParticleChange & G4ParticleChange::operator=(const G4ParticleChange &right)
{
   if (verboseLevel>1) {
    G4cout << "G4ParticleChange:: assignment operator is called " << G4endl;
   }
   if (this != &right)
   {
      theListOfSecondaries = right.theListOfSecondaries;
      theSizeOftheListOfSecondaries = right.theSizeOftheListOfSecondaries;
      theNumberOfSecondaries = right.theNumberOfSecondaries;
      theStatusChange = right.theStatusChange;
      theCurrentTrack = right.theCurrentTrack;

      theMomentumDirectionChange = right.theMomentumDirectionChange;
      thePolarizationChange = right.thePolarizationChange;
      thePositionChange = right.thePositionChange;
      theTimeChange = right.theTimeChange;
      theEnergyChange = right.theEnergyChange;
      theMassChange = right.theMassChange;
      theChargeChange = right.theChargeChange;
      theMagneticMomentChange = right.theMagneticMomentChange;
      // theWeightChange = right.theWeightChange;

      theTrueStepLength = right.theTrueStepLength;
      theLocalEnergyDeposit = right.theLocalEnergyDeposit;
      theSteppingControlFlag = right.theSteppingControlFlag;
   }
   return *this;
}

G4bool G4ParticleChange::operator==(const G4ParticleChange &right) const
{
   return ((G4VParticleChange *)this == (G4VParticleChange *) &right);
}

G4bool G4ParticleChange::operator!=(const G4ParticleChange &right) const
{
   return ((G4VParticleChange *)this != (G4VParticleChange *) &right);
}


//----------------------------------------------------------------
// methods for handling secondaries 
//

void G4ParticleChange::AddSecondary(G4DynamicParticle* aParticle, 
				    G4bool   IsGoodForTracking    )
{
  //  create track
  G4Track* aTrack = new G4Track(aParticle, theTimeChange, thePositionChange);

  // set IsGoodGorTrackingFlag
  if (IsGoodForTracking) aTrack->SetGoodForTrackingFlag();

  //   Touchable handle is copied to keep the pointer
  aTrack->SetTouchableHandle(theCurrentTrack->GetTouchableHandle());
 
 //  add a secondary
  G4VParticleChange::AddSecondary(aTrack);
}

void G4ParticleChange::AddSecondary(G4DynamicParticle* aParticle, 
				    G4ThreeVector      newPosition,
				    G4bool   IsGoodForTracking    )
{
  //  create track
  G4Track*  aTrack = new G4Track(aParticle, theTimeChange, newPosition);

  // set IsGoodGorTrackingFlag
  if (IsGoodForTracking) aTrack->SetGoodForTrackingFlag();

  //   Touchable is a temporary object, so you cannot keep the pointer
  aTrack->SetTouchableHandle((G4VTouchable*)0);

  //  add a secondary
  G4VParticleChange::AddSecondary(aTrack);
}

void G4ParticleChange::AddSecondary(G4DynamicParticle* aParticle, 
				    G4double           newTime,
				    G4bool   IsGoodForTracking    )
{
  //  create track
  G4Track*  aTrack = new G4Track(aParticle, newTime, thePositionChange); 

  // set IsGoodGorTrackingFlag
  if (IsGoodForTracking) aTrack->SetGoodForTrackingFlag();
 
  //   Touchable handle is copied to keep the pointer
  aTrack->SetTouchableHandle(theCurrentTrack->GetTouchableHandle());

  //  add a secondary
  G4VParticleChange::AddSecondary(aTrack);
}

void G4ParticleChange::AddSecondary(G4Track* aTrack)
{
  //  add a secondary
  G4VParticleChange::AddSecondary(aTrack);
}

//----------------------------------------------------------------
// functions for Initialization
//

void G4ParticleChange::Initialize(const G4Track& track)
{
  // use base class's method at first
  G4VParticleChange::Initialize(track);
  theCurrentTrack= &track;

  // set Energy/Momentum etc. equal to those of the parent particle
  const G4DynamicParticle*  pParticle = track.GetDynamicParticle();
  theEnergyChange          = pParticle->GetKineticEnergy();
  theMomentumDirectionChange        = pParticle->GetMomentumDirection();
  thePolarizationChange    = pParticle->GetPolarization();
  theProperTimeChange      = pParticle->GetProperTime();

  // Set mass/charge/MagneticMoment  of DynamicParticle
  theMassChange = pParticle->GetMass();
  theChargeChange = pParticle->GetCharge();
  theMagneticMomentChange = pParticle->GetMagneticMoment();

  // set Position/Time etc. equal to those of the parent track
  thePositionChange      = track.GetPosition();
  theTimeChange          = track.GetGlobalTime();

  // theWeightChange        = track.GetWeight();
}

//----------------------------------------------------------------
// methods for updating G4Step 
//

G4Step* G4ParticleChange::UpdateStepForAlongStep(G4Step* pStep)
{
  // A physics process always calculates the final state of the
  // particle relative to the initial state at the beginning
  // of the Step, i.e., based on information of G4Track (or
  // equivalently the PreStepPoint). 
  // So, the differences (delta) between these two states have to be
  // calculated and be accumulated in PostStepPoint. 
  
  // Take note that the return type of GetMomentumDirectionChange is a
  // pointer to G4ParticleMometum. Also it is a normalized 
  // momentum vector.

  G4StepPoint* pPreStepPoint  = pStep->GetPreStepPoint(); 
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint(); 
  G4double     mass = theMassChange;

  // Set Mass/Charge/MagneticMoment 
  pPostStepPoint->SetMass(theMassChange);
  pPostStepPoint->SetCharge(theChargeChange);  
  pPostStepPoint->SetMagneticMoment(theMagneticMomentChange);  
 
  // calculate new kinetic energy
  G4double energy = pPostStepPoint->GetKineticEnergy() 
                    + (theEnergyChange - pPreStepPoint->GetKineticEnergy()); 

  // update kinetic energy and momentum direction
  if (energy > 0.0) {
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
  } else {
    // stop case
    //pPostStepPoint->SetMomentumDirection(G4ThreeVector(1., 0., 0.));
    pPostStepPoint->SetKineticEnergy(0.0);
  }

  // update polarization
  pPostStepPoint->AddPolarization( thePolarizationChange
				   - pPreStepPoint->GetPolarization());
      
  // update position and time
  pPostStepPoint->AddPosition( thePositionChange 
			       - pPreStepPoint->GetPosition() );
  pPostStepPoint->AddGlobalTime( theTimeChange
				 - pPreStepPoint->GetGlobalTime());
  pPostStepPoint->AddLocalTime( theTimeChange 
				 - pPreStepPoint->GetGlobalTime()); 
  pPostStepPoint->AddProperTime( theProperTimeChange 
				 - pPreStepPoint->GetProperTime());

  if (isParentWeightProposed) {
    // update weight
    G4double newWeight= theParentWeight/(pPreStepPoint->GetWeight())
      * (pPostStepPoint->GetWeight());
    if (isParentWeightSetByProcess) pPostStepPoint->SetWeight( newWeight );
    
    if (!fSetSecondaryWeightByProcess) {    
      // Set weight of secondary tracks
      for (G4int index= 0; index<theNumberOfSecondaries; index++){
        if ( (*theListOfSecondaries)[index] ) {
          ((*theListOfSecondaries)[index])->SetWeight(newWeight); ;
        }
      }
    }
  }

#ifdef G4VERBOSE
  G4Track*     aTrack  = pStep->GetTrack();
  if (debugFlag) CheckIt(*aTrack);
#endif

  //  Update the G4Step specific attributes 
  return UpdateStepInfo(pStep);
}

G4Step* G4ParticleChange::UpdateStepForPostStep(G4Step* pStep)
{ 
  // A physics process always calculates the final state of the particle

  // Take note that the return type of GetMomentumChange is a
  // pointer to G4ParticleMometum. Also it is a normalized 
  // momentum vector.

  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint(); 
  G4Track*     aTrack  = pStep->GetTrack();

  // Set Mass/Charge
  pPostStepPoint->SetMass(theMassChange);
  pPostStepPoint->SetCharge(theChargeChange);  
  pPostStepPoint->SetMagneticMoment(theMagneticMomentChange);  
 
  // update kinetic energy and momentum direction
  pPostStepPoint->SetMomentumDirection(theMomentumDirectionChange);
  pPostStepPoint->SetKineticEnergy( theEnergyChange );

   // update polarization
  pPostStepPoint->SetPolarization( thePolarizationChange );
      
  // update position and time
  pPostStepPoint->SetPosition( thePositionChange  );
  pPostStepPoint->SetGlobalTime( theTimeChange  );
  pPostStepPoint->AddLocalTime( theTimeChange 
				 - aTrack->GetGlobalTime());
  pPostStepPoint->SetProperTime( theProperTimeChange  );

  if (isParentWeightProposed) {
    // update weight
    if( isParentWeightSetByProcess) pPostStepPoint->SetWeight( theParentWeight );
    if (!fSetSecondaryWeightByProcess) {    
      // Set weight of secondary tracks
      for (G4int index= 0; index<theNumberOfSecondaries; index++){
        if ( (*theListOfSecondaries)[index] ) {
          ((*theListOfSecondaries)[index])->SetWeight(theParentWeight); ;
        }
      }
    }
  }

#ifdef G4VERBOSE
  if (debugFlag) CheckIt(*aTrack);
#endif

  //  Update the G4Step specific attributes 
  return UpdateStepInfo(pStep);
}


G4Step* G4ParticleChange::UpdateStepForAtRest(G4Step* pStep)
{ 
  // A physics process always calculates the final state of the particle

  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint(); 
  G4Track*     aTrack  = pStep->GetTrack();

  // Set Mass/Charge
  pPostStepPoint->SetMass(theMassChange);
  pPostStepPoint->SetCharge(theChargeChange);  
  pPostStepPoint->SetMagneticMoment(theMagneticMomentChange);  
 
  // update kinetic energy and momentum direction
  pPostStepPoint->SetMomentumDirection(theMomentumDirectionChange);
  pPostStepPoint->SetKineticEnergy( theEnergyChange );

  // update polarization
  pPostStepPoint->SetPolarization( thePolarizationChange );
      
  // update position and time
  pPostStepPoint->SetPosition( thePositionChange  );
  pPostStepPoint->SetGlobalTime( theTimeChange  );
  pPostStepPoint->AddLocalTime( theTimeChange 
				 - aTrack->GetGlobalTime());
  pPostStepPoint->SetProperTime( theProperTimeChange  );

  if (isParentWeightProposed ) {
    // update weight 
    if (isParentWeightSetByProcess) pPostStepPoint->SetWeight( theParentWeight );
    if (!fSetSecondaryWeightByProcess) {    
      // Set weight of secondary tracks
      for (G4int index= 0; index<theNumberOfSecondaries; index++){
        if ( (*theListOfSecondaries)[index] ) {
          ((*theListOfSecondaries)[index])->SetWeight(theParentWeight); ;
        }
      }
    }
  }

#ifdef G4VERBOSE
  if (debugFlag) CheckIt(*aTrack);
#endif

  //  Update the G4Step specific attributes 
  return UpdateStepInfo(pStep);
}

//----------------------------------------------------------------
// methods for printing messages  
//

void G4ParticleChange::DumpInfo() const
{
// use base-class DumpInfo
  G4VParticleChange::DumpInfo();

  G4int oldprc = G4cout.precision(3);

  G4cout << "        Mass (GeV)   : " 
       << std::setw(20) << theMassChange/GeV
       << G4endl; 
  G4cout << "        Charge (eplus)   : " 
       << std::setw(20) << theChargeChange/eplus
       << G4endl; 
  G4cout << "        MagneticMoment   : " 
         << std::setw(20) << theMagneticMomentChange << G4endl;
  G4cout << "                :  = " << std::setw(20) 
         << theMagneticMomentChange*2.*theMassChange/c_squared/eplus/hbar_Planck 
         <<  "*[e hbar]/[2 m]" 
         << G4endl; 
  G4cout << "        Position - x (mm)   : " 
       << std::setw(20) << thePositionChange.x()/mm
       << G4endl; 
  G4cout << "        Position - y (mm)   : " 
       << std::setw(20) << thePositionChange.y()/mm
       << G4endl; 
  G4cout << "        Position - z (mm)   : " 
       << std::setw(20) << thePositionChange.z()/mm
       << G4endl;
  G4cout << "        Time (ns)           : " 
       << std::setw(20) << theTimeChange/ns
       << G4endl;
  G4cout << "        Proper Time (ns)    : " 
       << std::setw(20) << theProperTimeChange/ns
       << G4endl;
  G4cout << "        Momentum Direct - x : " 
       << std::setw(20) << theMomentumDirectionChange.x()
       << G4endl;
  G4cout << "        Momentum Direct - y : " 
       << std::setw(20) << theMomentumDirectionChange.y()
       << G4endl;
  G4cout << "        Momentum Direct - z : " 
       << std::setw(20) << theMomentumDirectionChange.z()
       << G4endl;
  G4cout << "        Kinetic Energy (MeV): " 
       << std::setw(20) << theEnergyChange/MeV
       << G4endl;
  G4cout << "        Polarization - x    : " 
       << std::setw(20) << thePolarizationChange.x()
       << G4endl;
  G4cout << "        Polarization - y    : " 
       << std::setw(20) << thePolarizationChange.y()
       << G4endl;
  G4cout << "        Polarization - z    : " 
       << std::setw(20) <<  thePolarizationChange.z()
       << G4endl;
  G4cout.precision(oldprc);
}

G4bool G4ParticleChange::CheckIt(const G4Track& aTrack)
{
  G4bool    exitWithError = false;
  G4double  accuracy;

  // No check in case of "fStopAndKill" 
  if (GetTrackStatus() ==   fStopAndKill )  {
    return G4VParticleChange::CheckIt(aTrack);
  }

  // MomentumDirection should be unit vector
  G4bool itsOKforMomentum = true;  
  if ( theEnergyChange >0.) {
    accuracy = std::fabs(theMomentumDirectionChange.mag2()-1.0);
    if (accuracy > accuracyForWarning) {
#ifdef G4VERBOSE
      G4cout << "  G4ParticleChange::CheckIt  : ";
      G4cout << "the Momentum Change is not unit vector !!" << G4endl;
      G4cout << "  Difference:  " << accuracy << G4endl;
#endif
      itsOKforMomentum = false;
      if (accuracy > accuracyForException) exitWithError = true;
    }
  }

  // Both global and proper time should not go back
  G4bool itsOKforGlobalTime = true;  
  accuracy = (aTrack.GetGlobalTime()- theTimeChange)/ns;
  if (accuracy > accuracyForWarning) {
#ifdef G4VERBOSE
    G4cout << "  G4ParticleChange::CheckIt    : ";
    G4cout << "the global time goes back  !!" << G4endl;
    G4cout << "  Difference:  " << accuracy  << "[ns] " <<G4endl;
#endif
    itsOKforGlobalTime = false;
    if (accuracy > accuracyForException) exitWithError = true;
  }

  G4bool itsOKforProperTime = true;
  accuracy = (aTrack.GetProperTime() - theProperTimeChange )/ns;
  if (accuracy > accuracyForWarning) {
#ifdef G4VERBOSE
    G4cout << "  G4ParticleChange::CheckIt    : ";
    G4cout << "the proper time goes back  !!" << G4endl;
    G4cout << "  Difference:  " << accuracy  << "[ns] " <<G4endl;
#endif
    itsOKforProperTime = false;
    if (accuracy > accuracyForException) exitWithError = true;
  }

  // Kinetic Energy should not be negative
  G4bool itsOKforEnergy = true;
  accuracy = -1.0*theEnergyChange/MeV;
  if (accuracy > accuracyForWarning) {
#ifdef G4VERBOSE
    G4cout << "  G4ParticleChange::CheckIt    : ";
    G4cout << "the kinetic energy is negative  !!" << G4endl;
    G4cout << "  Difference:  " << accuracy  << "[MeV] " <<G4endl;
#endif
    itsOKforEnergy = false;
    if (accuracy > accuracyForException) exitWithError = true;
  }

  G4bool itsOK = itsOKforMomentum && itsOKforEnergy && itsOKforProperTime && itsOKforGlobalTime;
  // dump out information of this particle change
#ifdef G4VERBOSE
  if (!itsOK) { 
    G4cout << " G4ParticleChange::CheckIt " <<G4endl;
    DumpInfo();
  }
#endif

  // Exit with error
  if (exitWithError) {
    G4Exception("G4ParticleChange::CheckIt",
		"200",
		EventMustBeAborted,
		"momentum, energy, and/or time was illegal");
  }
  //correction
  if (!itsOKforMomentum) {
    G4double vmag = theMomentumDirectionChange.mag();
    theMomentumDirectionChange = (1./vmag)*theMomentumDirectionChange;
  }
  if (!itsOKforGlobalTime) {
    theTimeChange = aTrack.GetGlobalTime();
  }
  if (!itsOKforProperTime) {
    theProperTimeChange = aTrack.GetProperTime();
  }
  if (!itsOKforEnergy) {
    theEnergyChange = 0.0;
  }

  itsOK = (itsOK) && G4VParticleChange::CheckIt(aTrack);
  return itsOK;
}



