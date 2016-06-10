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
// $Id: G4ParticleChange.cc 92776 2015-09-16 06:57:55Z gcosmo $
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
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4TrackFastVector.hh"
#include "G4DynamicParticle.hh"
#include "G4ExceptionSeverity.hh"

G4ParticleChange::G4ParticleChange()
  : G4VParticleChange(),  
    theMomentumDirectionChange(),
    thePolarizationChange(),
    theEnergyChange(0.), 
    theVelocityChange(0.), isVelocityChanged(false),
    thePositionChange(),
    theGlobalTime0(0.),    theLocalTime0(0.),
    theTimeChange(0.),     theProperTimeChange(0.), 
    theMassChange(0.), theChargeChange(0.),
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
G4ParticleChange::G4ParticleChange(const G4ParticleChange &right)
  : G4VParticleChange(right)
{
   if (verboseLevel>1) {
    G4cout << "G4ParticleChange::  copy constructor is called " << G4endl;
   }
   theCurrentTrack = right.theCurrentTrack;

   theMomentumDirectionChange = right.theMomentumDirectionChange;
   thePolarizationChange = right.thePolarizationChange;
   thePositionChange   = right.thePositionChange;
   theGlobalTime0      = right.theGlobalTime0;
   theLocalTime0       = right.theLocalTime0;
   theTimeChange       = right.theTimeChange;
   theProperTimeChange = right.theProperTimeChange;
   theEnergyChange     = right.theEnergyChange;
   theVelocityChange   = right.theVelocityChange;
   isVelocityChanged   = true;
   theMassChange       = right.theMassChange;
   theChargeChange     = right.theChargeChange;
   theMagneticMomentChange = right.theMagneticMomentChange;
}

// assignemnt operator
G4ParticleChange & G4ParticleChange::operator=(const G4ParticleChange &right)
{
#ifdef G4VERBOSE
   if (verboseLevel>1) {
    G4cout << "G4ParticleChange:: assignment operator is called " << G4endl;
   }
#endif
   if (this != &right){
     if (theNumberOfSecondaries>0) {
#ifdef G4VERBOSE
       if (verboseLevel>0) {
	 G4cout << "G4ParticleChange: assignment operator Warning  ";
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
     theCurrentTrack = right.theCurrentTrack;
     
     theMomentumDirectionChange = right.theMomentumDirectionChange;
     thePolarizationChange = right.thePolarizationChange;
     thePositionChange = right.thePositionChange;
     theGlobalTime0      = right.theGlobalTime0;
     theLocalTime0       = right.theLocalTime0;
     theTimeChange       = right.theTimeChange;
     theProperTimeChange = right.theProperTimeChange;
     theEnergyChange     = right.theEnergyChange;
     theVelocityChange   = right.theVelocityChange;
     isVelocityChanged   = true;
     theMassChange       = right.theMassChange;
     theChargeChange     = right.theChargeChange;
     theMagneticMomentChange = right.theMagneticMomentChange;
     
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
  G4Track* aTrack = new G4Track(aParticle, GetGlobalTime(), thePositionChange);

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
  G4Track*  aTrack = new G4Track(aParticle, GetGlobalTime(), newPosition);

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
  theEnergyChange            = pParticle->GetKineticEnergy();
  theVelocityChange          = track.GetVelocity();
  isVelocityChanged          = false;
  theMomentumDirectionChange = pParticle->GetMomentumDirection();
  thePolarizationChange      = pParticle->GetPolarization();
  theProperTimeChange        = pParticle->GetProperTime();

  // Set mass/charge/MagneticMoment  of DynamicParticle
  theMassChange = pParticle->GetMass();
  theChargeChange = pParticle->GetCharge();
  theMagneticMomentChange = pParticle->GetMagneticMoment();

  // set Position  equal to those of the parent track
  thePositionChange      = track.GetPosition();

  // set TimeChange equal to local time of the parent track
  theTimeChange                = track.GetLocalTime();

  // set initial Local/Global time of the parent track
  theLocalTime0           = track.GetLocalTime();
  theGlobalTime0          = track.GetGlobalTime();

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
  G4Track* pTrack = pStep->GetTrack();
  G4double     mass = theMassChange;

  // Set Mass/Charge/MagneticMoment 
  pPostStepPoint->SetMass(theMassChange);
  pPostStepPoint->SetCharge(theChargeChange);  
  pPostStepPoint->SetMagneticMoment(theMagneticMomentChange);  
 
  // calculate new kinetic energy
  G4double preEnergy = pPreStepPoint->GetKineticEnergy();
  G4double energy = pPostStepPoint->GetKineticEnergy() 
                    + (theEnergyChange - preEnergy); 

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
  // calculate velocity
  if (!isVelocityChanged) {
    if(energy > 0.0) {
      pTrack->SetKineticEnergy(energy);
      theVelocityChange = pTrack->CalculateVelocity();
      pTrack->SetKineticEnergy(preEnergy);
    } else if(theMassChange > 0.0) {
      theVelocityChange = 0.0;
    }
  }
  pPostStepPoint->SetVelocity(theVelocityChange);

  // update polarization
  pPostStepPoint->AddPolarization( thePolarizationChange
				   - pPreStepPoint->GetPolarization());
      
  // update position and time
  pPostStepPoint->AddPosition( thePositionChange 
			       - pPreStepPoint->GetPosition() );
  pPostStepPoint->AddGlobalTime(theTimeChange - theLocalTime0);
  pPostStepPoint->AddLocalTime( theTimeChange - theLocalTime0 );	       
  pPostStepPoint->AddProperTime( theProperTimeChange 
				 - pPreStepPoint->GetProperTime());

  if (isParentWeightProposed ){
    pPostStepPoint->SetWeight( theParentWeight );
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
  G4Track* pTrack = pStep->GetTrack();

  // Set Mass/Charge
  pPostStepPoint->SetMass(theMassChange);
  pPostStepPoint->SetCharge(theChargeChange);  
  pPostStepPoint->SetMagneticMoment(theMagneticMomentChange);  
 
  // update kinetic energy and momentum direction
  pPostStepPoint->SetMomentumDirection(theMomentumDirectionChange);
  pPostStepPoint->SetKineticEnergy( theEnergyChange );

  // calculate velocity
  pTrack->SetKineticEnergy( theEnergyChange );
  if (!isVelocityChanged) {
    if(theEnergyChange > 0.0) {
      theVelocityChange = pTrack->CalculateVelocity();
    } else if(theMassChange > 0.0) {
      theVelocityChange = 0.0;
    }
  }
  pPostStepPoint->SetVelocity(theVelocityChange);
 
   // update polarization
  pPostStepPoint->SetPolarization( thePolarizationChange );
      
  // update position and time
  pPostStepPoint->SetPosition( thePositionChange  );
  pPostStepPoint->AddGlobalTime(theTimeChange - theLocalTime0);
  pPostStepPoint->SetLocalTime( theTimeChange );	       
  pPostStepPoint->SetProperTime( theProperTimeChange  );

  if (isParentWeightProposed ){
    pPostStepPoint->SetWeight( theParentWeight );
  }

#ifdef G4VERBOSE
  G4Track*     aTrack  = pStep->GetTrack();
  if (debugFlag) CheckIt(*aTrack);
#endif

  //  Update the G4Step specific attributes 
  return UpdateStepInfo(pStep);
}


G4Step* G4ParticleChange::UpdateStepForAtRest(G4Step* pStep)
{ 
  // A physics process always calculates the final state of the particle

  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint(); 

  // Set Mass/Charge
  pPostStepPoint->SetMass(theMassChange);
  pPostStepPoint->SetCharge(theChargeChange);  
  pPostStepPoint->SetMagneticMoment(theMagneticMomentChange);  
 
  // update kinetic energy and momentum direction
  pPostStepPoint->SetMomentumDirection(theMomentumDirectionChange);
  pPostStepPoint->SetKineticEnergy( theEnergyChange );
  if (!isVelocityChanged) theVelocityChange = pStep->GetTrack()->CalculateVelocity();
  pPostStepPoint->SetVelocity(theVelocityChange);

  // update polarization
  pPostStepPoint->SetPolarization( thePolarizationChange );
      
  // update position and time
  pPostStepPoint->SetPosition( thePositionChange  );
  pPostStepPoint->AddGlobalTime(theTimeChange - theLocalTime0);
  pPostStepPoint->SetLocalTime( theTimeChange );	       
  pPostStepPoint->SetProperTime( theProperTimeChange  );

  if (isParentWeightProposed ){
    pPostStepPoint->SetWeight( theParentWeight );
  }

#ifdef G4VERBOSE
  G4Track*     aTrack  = pStep->GetTrack();
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
  G4cout << "        Velocity  (/c): " 
       << std::setw(20) << theVelocityChange/c_light
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
  static G4ThreadLocal G4int nError = 0;
#ifdef G4VERBOSE
  const  G4int maxError = 30;
#endif

  // No check in case of "fStopAndKill" 
  if (GetTrackStatus() ==   fStopAndKill )  return G4VParticleChange::CheckIt(aTrack);

  // MomentumDirection should be unit vector
  G4bool itsOKforMomentum = true;  
  if ( theEnergyChange >0.) {
    accuracy = std::fabs(theMomentumDirectionChange.mag2()-1.0);
    if (accuracy > accuracyForWarning) {
      itsOKforMomentum = false;
      nError += 1;
      exitWithError = exitWithError || (accuracy > accuracyForException);
#ifdef G4VERBOSE
      if (nError < maxError) {
	G4cout << "  G4ParticleChange::CheckIt  : ";
	G4cout << "the Momentum Change is not unit vector !!" 
	       << "  Difference:  " << accuracy << G4endl;
	G4cout << aTrack.GetDefinition()->GetParticleName()
	       << " E=" << aTrack.GetKineticEnergy()/MeV
	       << " pos=" << aTrack.GetPosition().x()/m
	       << ", " << aTrack.GetPosition().y()/m
	       << ", " << aTrack.GetPosition().z()/m
	       <<G4endl;
      }
#endif
    }
  }

  // Both global and proper time should not go back
  G4bool itsOKforGlobalTime = true;  
  accuracy = (aTrack.GetLocalTime()- theTimeChange)/ns;
  if (accuracy > accuracyForWarning) {
    itsOKforGlobalTime = false;
    nError += 1;
    exitWithError = exitWithError || (accuracy > accuracyForException);
#ifdef G4VERBOSE
    if (nError < maxError) {
      G4cout << "  G4ParticleChange::CheckIt    : ";
      G4cout << "the local time goes back  !!" 
	     << "  Difference:  " << accuracy  << "[ns] " <<G4endl;
      G4cout << aTrack.GetDefinition()->GetParticleName()
	     << " E=" << aTrack.GetKineticEnergy()/MeV
	     << " pos=" << aTrack.GetPosition().x()/m
	     << ", " << aTrack.GetPosition().y()/m
	     << ", " << aTrack.GetPosition().z()/m
	     << " global time=" << aTrack.GetGlobalTime()/ns
	     << " local time=" << aTrack.GetLocalTime()/ns
	     << " proper time=" << aTrack.GetProperTime()/ns
	     << G4endl;
    }
#endif
  }

  G4bool itsOKforProperTime = true;
  accuracy = (aTrack.GetProperTime() - theProperTimeChange )/ns;
  if (accuracy > accuracyForWarning) {
    itsOKforProperTime = false;
    nError += 1;
    exitWithError = exitWithError ||  (accuracy > accuracyForException);
#ifdef G4VERBOSE
    if (nError < maxError) {
      G4cout << "  G4ParticleChange::CheckIt    : ";
      G4cout << "the proper time goes back  !!" 
	     << "  Difference:  " << accuracy  << "[ns] " <<G4endl;
      G4cout << aTrack.GetDefinition()->GetParticleName()
	     << " E=" << aTrack.GetKineticEnergy()/MeV
	     << " pos=" << aTrack.GetPosition().x()/m
	     << ", " << aTrack.GetPosition().y()/m
	     << ", " << aTrack.GetPosition().z()/m
	     << " global time=" << aTrack.GetGlobalTime()/ns
	     << " local time=" << aTrack.GetLocalTime()/ns
	     << " proper time=" << aTrack.GetProperTime()/ns
	     <<G4endl;
    }
#endif
  }

  // Kinetic Energy should not be negative
  G4bool itsOKforEnergy = true;
  accuracy = -1.0*theEnergyChange/MeV;
  if (accuracy > accuracyForWarning) {
    itsOKforEnergy = false;
    nError += 1;
    exitWithError = exitWithError ||   (accuracy > accuracyForException);
#ifdef G4VERBOSE
    if (nError < maxError) {
      G4cout << "  G4ParticleChange::CheckIt    : ";
      G4cout << "the kinetic energy is negative  !!" 
	     << "  Difference:  " << accuracy  << "[MeV] " <<G4endl;
      G4cout << aTrack.GetDefinition()->GetParticleName()
	     << " E=" << aTrack.GetKineticEnergy()/MeV
	     << " pos=" << aTrack.GetPosition().x()/m
	     << ", " << aTrack.GetPosition().y()/m
	     << ", " << aTrack.GetPosition().z()/m
	     <<G4endl;
    }
#endif
  }

  // Velocity  should not be less than c_light
  G4bool itsOKforVelocity = true;
  if (theVelocityChange < 0.) {
    itsOKforVelocity = false;
    nError += 1;
    exitWithError = true;
#ifdef G4VERBOSE
    if (nError < maxError) {
      G4cout << "  G4ParticleChange::CheckIt    : ";
      G4cout << "the velocity is negative  !!" 
	     << "  Velocity:  " << theVelocityChange/c_light  <<G4endl;
      G4cout << aTrack.GetDefinition()->GetParticleName()
	     << " E=" << aTrack.GetKineticEnergy()/MeV
	     << " pos=" << aTrack.GetPosition().x()/m
	     << ", " << aTrack.GetPosition().y()/m
	     << ", " << aTrack.GetPosition().z()/m
	     <<G4endl;
    }
#endif
  }

  accuracy = theVelocityChange/c_light - 1.0;
  if (accuracy > accuracyForWarning) {
    itsOKforVelocity = false;
    nError += 1;
    exitWithError = exitWithError ||  (accuracy > accuracyForException);
#ifdef G4VERBOSE
    if (nError < maxError) {
      G4cout << "  G4ParticleChange::CheckIt    : ";
      G4cout << "the velocity is greater than c_light  !!" << G4endl;
      G4cout << "  Velocity:  " << theVelocityChange/c_light  <<G4endl;
      G4cout << aTrack.GetDefinition()->GetParticleName()
	     << " E=" << aTrack.GetKineticEnergy()/MeV
	     << " pos=" << aTrack.GetPosition().x()/m
	     << ", " << aTrack.GetPosition().y()/m
	     << ", " << aTrack.GetPosition().z()/m
	     <<G4endl;
    }
#endif
  }

  G4bool itsOK = itsOKforMomentum && itsOKforEnergy && itsOKforVelocity && itsOKforProperTime && itsOKforGlobalTime;
  // dump out information of this particle change
#ifdef G4VERBOSE
  if (!itsOK) { 
    DumpInfo();
  }
#endif

  // Exit with error
  if (exitWithError) {
    G4Exception("G4ParticleChange::CheckIt",
		"TRACK003", EventMustBeAborted,
		"momentum, energy, and/or time was illegal");
  }
  //correction
  if (!itsOKforMomentum) {
    G4double vmag = theMomentumDirectionChange.mag();
    theMomentumDirectionChange = (1./vmag)*theMomentumDirectionChange;
  }
  if (!itsOKforGlobalTime) {
    theTimeChange = aTrack.GetLocalTime();
  }
  if (!itsOKforProperTime) {
    theProperTimeChange = aTrack.GetProperTime();
  }
  if (!itsOKforEnergy) {
    theEnergyChange = 0.0;
  }
  if (!itsOKforVelocity) {
    theVelocityChange = c_light;
  }

  itsOK = (itsOK) && G4VParticleChange::CheckIt(aTrack);
  return itsOK;
}
