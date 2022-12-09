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
// G4ParticleChange class implementation
//
// Author: Hisaya Kurashige, 23 March 1998  
// --------------------------------------------------------------------

#include "G4ParticleChange.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4DynamicParticle.hh"
#include "G4ExceptionSeverity.hh"

// --------------------------------------------------------------------
G4ParticleChange::G4ParticleChange()
{}

// --------------------------------------------------------------------
void G4ParticleChange::AddSecondary(G4DynamicParticle* aParticle,
                                    G4bool IsGoodForTracking)
{
  // create track
  G4Track* aTrack = new G4Track(aParticle, GetGlobalTime(), thePositionChange);

  // set IsGoodGorTrackingFlag
  if(IsGoodForTracking)
    aTrack->SetGoodForTrackingFlag();

  // touchable handle is copied to keep the pointer
  aTrack->SetTouchableHandle(theCurrentTrack->GetTouchableHandle());

  // add a secondary
  G4VParticleChange::AddSecondary(aTrack);
}

// --------------------------------------------------------------------
void G4ParticleChange::AddSecondary(G4DynamicParticle* aParticle,
                                    G4ThreeVector newPosition,
                                    G4bool IsGoodForTracking)
{
  // create track
  G4Track* aTrack = new G4Track(aParticle, GetGlobalTime(), newPosition);

  // set IsGoodGorTrackingFlag
  if(IsGoodForTracking)
    aTrack->SetGoodForTrackingFlag();

  // touchable is a temporary object, so you cannot keep the pointer
  aTrack->SetTouchableHandle(nullptr);

  // add a secondary
  G4VParticleChange::AddSecondary(aTrack);
}

// --------------------------------------------------------------------
void G4ParticleChange::AddSecondary(G4DynamicParticle* aParticle,
                                    G4double newTime, G4bool IsGoodForTracking)
{
  // create track
  G4Track* aTrack = new G4Track(aParticle, newTime, thePositionChange);

  // set IsGoodGorTrackingFlag
  if(IsGoodForTracking)
    aTrack->SetGoodForTrackingFlag();

  // touchable handle is copied to keep the pointer
  aTrack->SetTouchableHandle(theCurrentTrack->GetTouchableHandle());

  // add a secondary
  G4VParticleChange::AddSecondary(aTrack);
}

// --------------------------------------------------------------------

void G4ParticleChange::AddSecondary(G4Track* aTrack)
{
  // add a secondary
  G4VParticleChange::AddSecondary(aTrack);
}

// --------------------------------------------------------------------
void G4ParticleChange::Initialize(const G4Track& track)
{
  // use base class's method at first
  G4VParticleChange::Initialize(track);

  // set Energy/Momentum etc. equal to those of the parent particle
  const G4DynamicParticle* pParticle = track.GetDynamicParticle();
  theEnergyChange                    = pParticle->GetKineticEnergy();
  theVelocityChange                  = track.GetVelocity();
  isVelocityChanged                  = false;
  theMomentumDirectionChange         = pParticle->GetMomentumDirection();
  thePolarizationChange              = pParticle->GetPolarization();
  theProperTimeChange                = pParticle->GetProperTime();

  // set mass/charge/MagneticMoment of DynamicParticle
  theMassChange           = pParticle->GetMass();
  theChargeChange         = pParticle->GetCharge();
  theMagneticMomentChange = pParticle->GetMagneticMoment();

  // set Position equal to those of the parent track
  thePositionChange = track.GetPosition();

  // set TimeChange equal to local time of the parent track
  theTimeChange = theLocalTime0 = track.GetLocalTime();

  // set initial global time of the parent track
  theGlobalTime0 = track.GetGlobalTime();
}

// --------------------------------------------------------------------
G4Step* G4ParticleChange::UpdateStepForAlongStep(G4Step* pStep)
{
  // A physics process always calculates the final state of the
  // particle relative to the initial state at the beginning
  // of the Step, i.e., based on information of G4Track (or
  // equivalently the PreStepPoint).
  // So, the differences (delta) between these two states have to be
  // calculated and be accumulated in PostStepPoint

  // Take note that the return type of GetMomentumDirectionChange() is a
  // pointer to G4ParticleMometum. Also it is a normalized
  // momentum vector

  const G4StepPoint* pPreStepPoint = pStep->GetPreStepPoint();
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();

  // set Mass/Charge/MagneticMoment
  pPostStepPoint->SetMass(theMassChange);
  pPostStepPoint->SetCharge(theChargeChange);
  pPostStepPoint->SetMagneticMoment(theMagneticMomentChange);

  // calculate new kinetic energy
  G4double preEnergy = pPreStepPoint->GetKineticEnergy();
  G4double energy =
    pPostStepPoint->GetKineticEnergy() + (theEnergyChange - preEnergy);

  // update kinetic energy and momentum direction
  if(energy > 0.0)
  {
    // calculate new momentum
    G4ThreeVector pMomentum = pPostStepPoint->GetMomentum()
      + (CalcMomentum(theEnergyChange, theMomentumDirectionChange,
		      theMassChange) - pPreStepPoint->GetMomentum());
    G4double tMomentum2 = pMomentum.mag2();
    G4ThreeVector direction(1.0, 0.0, 0.0);
    if(tMomentum2 > 0.)
    {
      direction = pMomentum / std::sqrt(tMomentum2);
    }
    pPostStepPoint->SetMomentumDirection(direction);
    pPostStepPoint->SetKineticEnergy(energy);

    // if velocity is not set it should be recomputed
    //  
    if(!isVelocityChanged)
    {
      if (theMassChange > 0.0)
      {
	theVelocityChange = CLHEP::c_light *
	  std::sqrt(energy*(energy + 2*theMassChange))/(energy + theMassChange);
      }
      else 
      {
	// zero mass particle
	theVelocityChange = CLHEP::c_light;
	// optical photon case
	if(theCurrentTrack->GetParticleDefinition()->GetPDGEncoding() == -22)
	{
          G4Track* pTrack = pStep->GetTrack();
          G4double e = pTrack->GetKineticEnergy();
	  pTrack->SetKineticEnergy(energy);
          theVelocityChange = pTrack->CalculateVelocityForOpticalPhoton();
	  pTrack->SetKineticEnergy(e);
	}
      }
    }
    pPostStepPoint->SetVelocity(theVelocityChange);
  }
  else
  {
    // stop case
    pPostStepPoint->SetKineticEnergy(0.0);
    pPostStepPoint->SetVelocity(0.0);
  }

  // update polarization
  pPostStepPoint->AddPolarization(thePolarizationChange -
                                  pPreStepPoint->GetPolarization());

  // update position and time
  pPostStepPoint->AddPosition(thePositionChange - pPreStepPoint->GetPosition());
  pPostStepPoint->AddGlobalTime(theTimeChange - theLocalTime0);
  pPostStepPoint->AddLocalTime(theTimeChange - theLocalTime0);
  pPostStepPoint->AddProperTime(theProperTimeChange -
                                pPreStepPoint->GetProperTime());

  if(isParentWeightProposed)
  {
    pPostStepPoint->SetWeight(theParentWeight);
  }

#ifdef G4VERBOSE
  if(debugFlag) { CheckIt( *theCurrentTrack ); }
#endif

  // update the G4Step specific attributes
  return UpdateStepInfo(pStep);
}

// --------------------------------------------------------------------
G4Step* G4ParticleChange::UpdateStepForPostStep(G4Step* pStep)
{
  // A physics process always calculates the final state of the particle

  // Take note that the return type of GetMomentumChange() is a
  // pointer to G4ParticleMometum. Also it is a normalized
  // momentum vector

  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();
  G4Track* pTrack             = pStep->GetTrack();

  // set Mass/Charge
  pPostStepPoint->SetMass(theMassChange);
  pPostStepPoint->SetCharge(theChargeChange);
  pPostStepPoint->SetMagneticMoment(theMagneticMomentChange);

  // update kinetic energy and momentum direction
  pPostStepPoint->SetMomentumDirection(theMomentumDirectionChange);

  // calculate velocity
  if(theEnergyChange > 0.0)
  {
    pPostStepPoint->SetKineticEnergy(theEnergyChange);
    pTrack->SetKineticEnergy(theEnergyChange);
    if(!isVelocityChanged)
    {
      theVelocityChange = pTrack->CalculateVelocity();
    }
    pPostStepPoint->SetVelocity(theVelocityChange);
  }
  else 
  {
    pPostStepPoint->SetKineticEnergy(0.0);
    pPostStepPoint->SetVelocity(0.0);
  }

  // update polarization
  pPostStepPoint->SetPolarization(thePolarizationChange);

  // update position and time
  pPostStepPoint->SetPosition(thePositionChange);
  pPostStepPoint->AddGlobalTime(theTimeChange - theLocalTime0);
  pPostStepPoint->SetLocalTime(theTimeChange);
  pPostStepPoint->SetProperTime(theProperTimeChange);

  if(isParentWeightProposed)
  {
    pPostStepPoint->SetWeight(theParentWeight);
  }

#ifdef G4VERBOSE
  if(debugFlag) { CheckIt( *theCurrentTrack ); }
#endif

  // update the G4Step specific attributes
  return UpdateStepInfo(pStep);
}

// --------------------------------------------------------------------
G4Step* G4ParticleChange::UpdateStepForAtRest(G4Step* pStep)
{
  // A physics process always calculates the final state of the particle

  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();

  // set Mass/Charge
  pPostStepPoint->SetMass(theMassChange);
  pPostStepPoint->SetCharge(theChargeChange);
  pPostStepPoint->SetMagneticMoment(theMagneticMomentChange);

  // update kinetic energy and momentum direction
  pPostStepPoint->SetMomentumDirection(theMomentumDirectionChange);
  pPostStepPoint->SetKineticEnergy(theEnergyChange);
  if(!isVelocityChanged)
    theVelocityChange = theCurrentTrack->CalculateVelocity();
  pPostStepPoint->SetVelocity(theVelocityChange);

  // update polarization
  pPostStepPoint->SetPolarization(thePolarizationChange);

  // update position and time
  pPostStepPoint->SetPosition(thePositionChange);
  pPostStepPoint->AddGlobalTime(theTimeChange - theLocalTime0);
  pPostStepPoint->SetLocalTime(theTimeChange);
  pPostStepPoint->SetProperTime(theProperTimeChange);

  if(isParentWeightProposed)
  {
    pPostStepPoint->SetWeight(theParentWeight);
  }

#ifdef G4VERBOSE
  if(debugFlag) { CheckIt( *theCurrentTrack ); }
#endif

  // update the G4Step specific attributes
  return UpdateStepInfo(pStep);
}

// --------------------------------------------------------------------
void G4ParticleChange::DumpInfo() const
{
  // use base-class DumpInfo
  G4VParticleChange::DumpInfo();

  G4long oldprc = G4cout.precision(8);

  G4cout << "        Mass (GeV)          : " << std::setw(20) << theMassChange / GeV
         << G4endl;
  G4cout << "        Charge (eplus)      : " << std::setw(20)
         << theChargeChange / eplus << G4endl;
  G4cout << "        MagneticMoment      : " << std::setw(20)
         << theMagneticMomentChange << G4endl;
  G4cout << "                         =  : " << std::setw(20)
         << theMagneticMomentChange
            * 2. * theMassChange / c_squared / eplus / hbar_Planck
         << "*[e hbar]/[2 m]" << G4endl;
  G4cout << "        Position - x (mm)   : " << std::setw(20)
         << thePositionChange.x() / mm << G4endl;
  G4cout << "        Position - y (mm)   : " << std::setw(20)
         << thePositionChange.y() / mm << G4endl;
  G4cout << "        Position - z (mm)   : " << std::setw(20)
         << thePositionChange.z() / mm << G4endl;
  G4cout << "        Time (ns)           : " << std::setw(20)
         << theTimeChange / ns << G4endl;
  G4cout << "        Proper Time (ns)    : " << std::setw(20)
         << theProperTimeChange / ns << G4endl;
  G4cout << "        Momentum Direct - x : " << std::setw(20)
         << theMomentumDirectionChange.x() << G4endl;
  G4cout << "        Momentum Direct - y : " << std::setw(20)
         << theMomentumDirectionChange.y() << G4endl;
  G4cout << "        Momentum Direct - z : " << std::setw(20)
         << theMomentumDirectionChange.z() << G4endl;
  G4cout << "        Kinetic Energy (MeV): " << std::setw(20)
         << theEnergyChange / MeV << G4endl;
  G4cout << "        Velocity  (/c)      : " << std::setw(20)
         << theVelocityChange / c_light << G4endl;
  G4cout << "        Polarization - x    : " << std::setw(20)
         << thePolarizationChange.x() << G4endl;
  G4cout << "        Polarization - y    : " << std::setw(20)
         << thePolarizationChange.y() << G4endl;
  G4cout << "        Polarization - z    : " << std::setw(20)
         << thePolarizationChange.z() << G4endl;
  G4cout.precision(oldprc);
}
