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
// G4ParticleChangeForTransport class implementation
//
// Author: Hisaya Kurashige, 10 May 1998
// --------------------------------------------------------------------

#include "G4ParticleChangeForTransport.hh"
#include "G4TouchableHandle.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4TrackFastVector.hh"
#include "G4DynamicParticle.hh"

// --------------------------------------------------------------------
G4ParticleChangeForTransport::G4ParticleChangeForTransport() 
{
  // Disable flag that is enabled in G4VParticleChange if G4VERBOSE.
  debugFlag = false;
}

// --------------------------------------------------------------------
G4Step* G4ParticleChangeForTransport::UpdateStepForAtRest(G4Step* pStep)
{
  // Update the G4Step specific attributes
  return UpdateStepInfo(pStep);
}

// --------------------------------------------------------------------
G4Step* G4ParticleChangeForTransport::UpdateStepForAlongStep(G4Step* pStep)
{
  // Smooth curved tajectory representation: let the Step know about
  // the auxiliary trajectory points (jacek 30/10/2002)
  pStep->SetPointerToVectorOfAuxiliaryPoints(fpVectorOfAuxiliaryPointsPointer);

  // Most of the code assumes that transportation is always the first process,
  // so the pre- and post-step point are still equal.
  G4StepPoint* pPreStepPoint  = pStep->GetPreStepPoint();
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();

  // update momentum direction and energy
  if(isMomentumChanged)
  {
    pPostStepPoint->SetMomentumDirection(theMomentumDirectionChange);
    pPostStepPoint->SetKineticEnergy(theEnergyChange);
  }
  if(isVelocityChanged)
    pPostStepPoint->SetVelocity(theVelocityChange);

  // update polarization
  pPostStepPoint->SetPolarization(thePolarizationChange);

  // update position and time
  pPostStepPoint->SetPosition(thePositionChange);
  pPostStepPoint->AddGlobalTime(theTimeChange - pPreStepPoint->GetLocalTime());
  pPostStepPoint->AddLocalTime(theTimeChange - pPreStepPoint->GetLocalTime());
  pPostStepPoint->SetProperTime(theProperTimeChange);

#ifdef G4VERBOSE
  if(debugFlag) { CheckIt(*theCurrentTrack); }
#endif

  // Update the G4Step specific attributes
  pStep->SetStepLength( theTrueStepLength );
  pStep->SetControlFlag(theSteppingControlFlag);

  return pStep;
}

// --------------------------------------------------------------------
G4Step* G4ParticleChangeForTransport::UpdateStepForPostStep(G4Step* pStep)
{
  // A physics process always calculates the final state of the particle

  // Change volume only if some kinetic energy remains
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();
  if(pPostStepPoint->GetKineticEnergy() > 0.0)
  {
    // update next touchable
    // (touchable can be changed only at PostStepDoIt)
    pPostStepPoint->SetTouchableHandle(theTouchableHandle);

    pPostStepPoint->SetMaterial(theMaterialChange);
    pPostStepPoint->SetMaterialCutsCouple(theMaterialCutsCoupleChange);
    pPostStepPoint->SetSensitiveDetector(theSensitiveDetectorChange);
  }
  if(this->GetFirstStepInVolume())
  {
    pStep->SetFirstStepFlag();
  }
  else
  {
    pStep->ClearFirstStepFlag();
  }
  if(this->GetLastStepInVolume())
  {
    pStep->SetLastStepFlag();
  }
  else
  {
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

// --------------------------------------------------------------------
void G4ParticleChangeForTransport::DumpInfo() const
{
  // use base-class DumpInfo
  G4ParticleChange::DumpInfo();
  G4cout << "        Touchable (pointer) : " << std::setw(20)
         << theTouchableHandle() << G4endl;
}
