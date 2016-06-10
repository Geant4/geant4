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
// $Id: G4StepPoint.hh 68795 2013-04-05 13:24:46Z gcosmo $
//
//
//---------------------------------------------------------------
//
// G4StepPoint.hh
//
// Class Description:
//   This class represents information associated with the
//   each end of a Step like the space/time data of the
//   particle.
//
// Contact:
//   Questions and comments to this code should be sent to
//     Hisaya Kurashige  
//
// ---------------------------------------------------------------
//   Added fpMaterial                       16 FEb. 2000  H.Kurahige
//   Added fpMaterialCutsCouple              8 Oct. 2002  H.Kurahige
//   Added fMagneticMoment                      Mar 2007  H.Kurashige

#ifndef G4StepPoint_h
#define G4StepPoint_h 1

#include <cmath>                     // Include from 'system'
#include <CLHEP/Units/PhysicalConstants.h>

#include "globals.hh"                // Include from 'global'
#include "G4Allocator.hh"            // Include from 'global'
#include "G4ThreeVector.hh"          // Include from 'geometry'
#include "G4VPhysicalVolume.hh"      // Include from 'geometry'
class G4VProcess;
#include "G4SteppingControl.hh"
#include "G4StepStatus.hh"           // Include from 'track'
#include "G4TouchableHandle.hh"      // Include from 'geometry'
#include "G4Material.hh"
#include "G4LogicalVolume.hh"

class G4MaterialCutsCouple;
class G4VSensitiveDetector;
/////////////////
class G4StepPoint
///////////////// 
{

//--------
   public:


// Constructor/Destructor
   G4StepPoint();

   ~G4StepPoint(){}

// Copy Counstructor and assignment operator
   G4StepPoint(const G4StepPoint& );
   G4StepPoint & operator=(const G4StepPoint &);   

//--------
   
  public: // with description 

// Get/Set functions
  const G4ThreeVector& GetPosition() const;
  void SetPosition(const G4ThreeVector& aValue);
  void AddPosition(const G4ThreeVector& aValue);

  G4double GetLocalTime() const;
  void SetLocalTime(const G4double aValue);
  void AddLocalTime(const G4double aValue);
      // Time since the track is created.
     
   G4double GetGlobalTime() const;
   void SetGlobalTime(const G4double aValue);
   void AddGlobalTime(const G4double aValue);
      // Time since the event in which the track belongs is created.
     
   G4double GetProperTime() const;
   void SetProperTime(const G4double aValue);
   void AddProperTime(const G4double aValue);
      // Proper time of the particle.
     
   const G4ThreeVector& GetMomentumDirection() const;
   void SetMomentumDirection(const G4ThreeVector& aValue);
   void AddMomentumDirection(const G4ThreeVector& aValue);
     // Direction of momentum  (should be an unit vector)
    
   G4ThreeVector GetMomentum() const;
      // Total momentum of the track

  
   G4double GetTotalEnergy() const;
     // Total energy of the track

   G4double GetKineticEnergy() const;
   void SetKineticEnergy(const G4double aValue);
   void AddKineticEnergy(const G4double aValue);
     // Kinetic Energy of the track

   G4double GetVelocity() const;
   void SetVelocity(G4double v);
   //

   G4double GetBeta() const;
    // Velocity of the track in unit of c(light velocity)

   G4double GetGamma() const;
     // Gamma factor (1/sqrt[1-beta*beta]) of the track    

   G4VPhysicalVolume* GetPhysicalVolume() const;

   const G4VTouchable* GetTouchable() const;
   const G4TouchableHandle& GetTouchableHandle() const;
   void  SetTouchableHandle(const G4TouchableHandle& apValue);

   G4Material* GetMaterial() const;
   void SetMaterial(G4Material*);

   const G4MaterialCutsCouple* GetMaterialCutsCouple() const;
   void SetMaterialCutsCouple(const G4MaterialCutsCouple*);

   G4VSensitiveDetector* GetSensitiveDetector() const;
   void SetSensitiveDetector(G4VSensitiveDetector*);

   G4double GetSafety() const;
   void SetSafety(const G4double aValue);

   const G4ThreeVector& GetPolarization() const;
   void SetPolarization(const G4ThreeVector& aValue);
   void AddPolarization(const G4ThreeVector& aValue);

   G4StepStatus GetStepStatus() const;
   void SetStepStatus(const G4StepStatus aValue);

   const G4VProcess* GetProcessDefinedStep() const;
     // If the pointer is 0, this means the Step is defined
     // by the user defined limit in the current volume.
   void SetProcessDefinedStep(const G4VProcess* aValue);

  
   G4double GetMass() const;
   void SetMass(G4double value);

   G4double GetCharge() const;
   void SetCharge(G4double value);

   G4double GetMagneticMoment() const;
   void SetMagneticMoment(G4double value);

   void SetWeight(G4double aValue);
   G4double GetWeight() const; 

//---------
   private:
//---------

// Member data
   G4ThreeVector fPosition;
   G4double fGlobalTime;         
      // Time since event is created
   G4double fLocalTime;          
      // Time since track is created
   G4double fProperTime;          
      // Time since track is created (in rest frame of particle)
   G4ThreeVector fMomentumDirection;
   G4double fKineticEnergy;
   G4double fVelocity; 
      //  Momentum,energy and velocity
   G4TouchableHandle fpTouchable;
      //  Touchable Handle  
   G4Material* fpMaterial;
      // Material of the volmue
   const G4MaterialCutsCouple* fpMaterialCutsCouple;
      // MaterialCutsCouple of the volmue
   G4VSensitiveDetector* fpSensitiveDetector;
   G4double fSafety;
   G4ThreeVector fPolarization;
   G4StepStatus fStepStatus;
      // DoIt type which defined the current Step.
   const G4VProcess* fpProcessDefinedStep;
      // Process which defined the current Step.
   G4double fMass;
      // Dynamical mass of the particle
   G4double fCharge;
      // Dynamical Charge of the particle
   G4double fMagneticMoment;
      // Dynamical MagneticMoment of the particle
   G4double fWeight;
      // Track Weight
};

#include "G4StepPoint.icc"

#endif
