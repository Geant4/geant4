// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4StepPoint.hh,v 1.7 2001-02-17 11:25:10 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
// ---------------------------------------------------------------
//   Added fpMaterial                       16 FEb. 2000  H.Kurahige

#ifndef G4StepPoint_h
#define G4StepPoint_h 1

#include "globals.hh"                // Include from 'global'
#include "G4Allocator.hh"            // Include from 'global'
#include "G4ThreeVector.hh"          // Include from 'geometry'
#include "G4VPhysicalVolume.hh"      // Include from 'geometry'
class G4VProcess;
#include "G4SteppingControl.hh"
#include "G4StepStatus.hh"           // Include from 'track'
#include "G4VTouchable.hh"           // Include from 'geometry'
#include "G4Material.hh"
#include "G4LogicalVolume.hh"

/////////////////
class G4StepPoint
///////////////// 
{

//--------
   public:


// Constructor/Destructor
   G4StepPoint();

   ~G4StepPoint(){};

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
  void  SetTouchable(const G4VTouchable* apValue);

  G4Material* GetMaterial() const;
  void SetMaterial(G4Material*);

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
  void SetProcessDefinedStep(G4VProcess* aValue);

  
  G4double GetMass() const;
  void SetMass(G4double value);

  G4double GetCharge() const;
  void SetCharge(G4double value);

  void SetWeight(G4double aValue);
  G4double GetWeight() const;

  public:
  // copy constructor and assignment operaor as protected
  G4StepPoint(const G4StepPoint &right);
  G4StepPoint & operator=(const G4StepPoint &right);

  

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
   const G4VTouchable* fpTouchable;
   G4Material* fpMaterial;
      // Material of the volmue
   G4double fSafety;
   G4ThreeVector fPolarization;
   G4StepStatus fStepStatus;
      // DoIt type which defined the current Step.
   G4VProcess* fpProcessDefinedStep;
      // Process which defined the current Step.
   G4double fMass;
      // Dynamical mass of the particle
   G4double fCharge;
      // Dynamical Charge of the particle
   G4double fWeight;
      // Track Weight
};

#include "G4StepPoint.icc"
#endif

