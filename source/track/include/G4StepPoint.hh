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
// G4StepPoint
//
// Class description:
//
// This class represents information associated with each end
// of a Step like the space/time data of the particle.

// Author: Hisaya Kurashige, 16 February 2000
// --------------------------------------------------------------------
#ifndef G4StepPoint_hh
#define G4StepPoint_hh 1

#include <cmath>  // Include from 'system'
#include <CLHEP/Units/PhysicalConstants.h>

#include "globals.hh"            // Include from 'global'
#include "G4Allocator.hh"        // Include from 'global'
#include "G4ThreeVector.hh"      // Include from 'geometry'
#include "G4VPhysicalVolume.hh"  // Include from 'geometry'
#include "G4SteppingControl.hh"
#include "G4StepStatus.hh"       // Include from 'track'
#include "G4TouchableHandle.hh"  // Include from 'geometry'
#include "G4Material.hh"
#include "G4LogicalVolume.hh"

class G4VProcess;
class G4MaterialCutsCouple;
class G4VSensitiveDetector;

class G4StepPoint
{
  public:

    G4StepPoint() = default;
    ~G4StepPoint()= default;
      // Constructor/Destructor

    G4StepPoint(const G4StepPoint&) = default;
    G4StepPoint& operator=(const G4StepPoint&);
      // Copy Constructor and assignment operator

    inline const G4ThreeVector& GetPosition() const;
    inline void SetPosition(const G4ThreeVector& aValue);
    inline void AddPosition(const G4ThreeVector& aValue);
      // Position

    inline G4double GetLocalTime() const;
    inline void SetLocalTime(const G4double aValue);
    inline void AddLocalTime(const G4double aValue);
      // Time since the track is created

    inline G4double GetGlobalTime() const;
    inline void SetGlobalTime(const G4double aValue);
    inline void AddGlobalTime(const G4double aValue);
      // Time since the event in which the track belongs is created

    inline G4double GetProperTime() const;
    inline void SetProperTime(const G4double aValue);
    inline void AddProperTime(const G4double aValue);
      // Proper time of the particle

    inline const G4ThreeVector& GetMomentumDirection() const;
    inline void SetMomentumDirection(const G4ThreeVector& aValue);
    inline void AddMomentumDirection(const G4ThreeVector& aValue);
      // Direction of momentum  (should be an unit vector)

    inline G4ThreeVector GetMomentum() const;
      // Total momentum of the track

    inline G4double GetTotalEnergy() const;
      // Total energy of the track

    inline G4double GetKineticEnergy() const;
    inline void SetKineticEnergy(const G4double aValue);
    inline void AddKineticEnergy(const G4double aValue);
      // Kinetic Energy of the track

    inline G4double GetVelocity() const;
    inline void SetVelocity(G4double v);
      // Velocity

    inline G4double GetBeta() const;
      // Velocity of the track in unit of c(light velocity)

    inline G4double GetGamma() const;
      // Gamma factor (1/sqrt[1-beta*beta]) of the track

    inline G4VPhysicalVolume* GetPhysicalVolume() const;

    inline const G4VTouchable* GetTouchable() const;
    inline const G4TouchableHandle& GetTouchableHandle() const;
    inline void SetTouchableHandle(const G4TouchableHandle& apValue);

    inline G4Material* GetMaterial() const;
    inline void SetMaterial(G4Material*);

    inline const G4MaterialCutsCouple* GetMaterialCutsCouple() const;
    inline void SetMaterialCutsCouple(const G4MaterialCutsCouple*);

    inline G4VSensitiveDetector* GetSensitiveDetector() const;
    inline void SetSensitiveDetector(G4VSensitiveDetector*);

    inline G4double GetSafety() const;
    inline void SetSafety(const G4double aValue);

    inline const G4ThreeVector& GetPolarization() const;
    inline void SetPolarization(const G4ThreeVector& aValue);
    inline void AddPolarization(const G4ThreeVector& aValue);

    inline G4StepStatus GetStepStatus() const;
    inline void SetStepStatus(const G4StepStatus aValue);

    inline const G4VProcess* GetProcessDefinedStep() const;
      // If the pointer is 0, this means the Step is defined
      // by the user defined limit in the current volume
    inline void SetProcessDefinedStep(const G4VProcess* aValue);

    inline G4double GetMass() const;
    inline void SetMass(G4double value);

    inline G4double GetCharge() const;
    inline void SetCharge(G4double value);

    inline G4double GetMagneticMoment() const;
    inline void SetMagneticMoment(G4double value);

    inline void SetWeight(G4double aValue);
    inline G4double GetWeight() const;

  private:

    G4ThreeVector fPosition;
    G4double fGlobalTime = 0.0;
      // Time since event is created
    G4double fLocalTime = 0.0;
      // Time since track is created
    G4double fProperTime = 0.0;
      // Time since track is created (in rest frame of particle)
    G4ThreeVector fMomentumDirection;
    G4double fKineticEnergy = 0.0;
    G4double fVelocity = 0.0;
      // Momentum,energy and velocity
    G4TouchableHandle fpTouchable;
      // Touchable Handle
    G4Material* fpMaterial = nullptr;
      // Material of the volmue
    const G4MaterialCutsCouple* fpMaterialCutsCouple = nullptr;
      // MaterialCutsCouple of the volmue
    G4VSensitiveDetector* fpSensitiveDetector = nullptr;
    G4double fSafety = 0.0;
    G4ThreeVector fPolarization;
    G4StepStatus fStepStatus = fUndefined;
      // DoIt type which defined the current Step.
    const G4VProcess* fpProcessDefinedStep = nullptr;
      // Process which defined the current Step.
    G4double fMass = 0.0;
      // Dynamical mass of the particle
    G4double fCharge = 0.0;
      // Dynamical Charge of the particle
    G4double fMagneticMoment = 0.0;
      // Dynamical MagneticMoment of the particle
    G4double fWeight = 0.0;
      // Track Weight
};

#include "G4StepPoint.icc"

#endif
