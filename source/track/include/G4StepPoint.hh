// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4StepPoint.hh,v 1.4 1999-11-07 16:32:02 kurasige Exp $
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
   ~G4StepPoint();

//--------
   public: // with description 

// Get/Set functions
   inline const G4ThreeVector& GetPosition() const
   { return fPosition; }
   inline void SetPosition(const G4ThreeVector& aValue)
   { fPosition = aValue; }
   inline void AddPosition(const G4ThreeVector& aValue)
   { fPosition += aValue; }
     // Position where the track locates
    
   inline G4double GetLocalTime() const
   { return fLocalTime; }
   inline void SetLocalTime(const G4double aValue)
   { fLocalTime = aValue; }
   inline void AddLocalTime(const G4double aValue)
   { fLocalTime += aValue; }
      // Time since the track is created.
     
   inline G4double GetGlobalTime() const
   { return fGlobalTime; }
   inline void SetGlobalTime(const G4double aValue)
   { fGlobalTime = aValue; }
   inline void AddGlobalTime(const G4double aValue)
   { fGlobalTime += aValue; }
      // Time since the event in which the track belongs is created.
     
   inline G4double GetProperTime() const
   { return fProperTime; }
   inline void SetProperTime(const G4double aValue)
   { fProperTime = aValue; }
   inline void AddProperTime(const G4double aValue)
   { fProperTime += aValue; }
      // Proper time of the particle.
     
   inline const G4ThreeVector& GetMomentumDirection() const
   { return fMomentumDirection; }
   inline void SetMomentumDirection(const G4ThreeVector& aValue)
   { fMomentumDirection = aValue;
   }
   inline void AddMomentumDirection(const G4ThreeVector& aValue)
   { fMomentumDirection += aValue;
   }
     // Direction of momentum  (should be an unit vector)
    
   inline G4ThreeVector GetMomentum() const
   { 
     G4double tMomentum = sqrt(fKineticEnergy*fKineticEnergy +
                               2*fKineticEnergy*fMass);
     return G4ThreeVector(fMomentumDirection.x()*tMomentum,
			  fMomentumDirection.y()*tMomentum,
			  fMomentumDirection.z()*tMomentum);
   }
     // Total momentum of the track

   inline G4double GetTotalEnergy() const
   { 
     return fKineticEnergy + fMass; 
   }
     // Total energy of the track

   inline G4double GetKineticEnergy() const
   { return fKineticEnergy; }
   inline void SetKineticEnergy(const G4double aValue)
   { fKineticEnergy = aValue; }
   inline void AddKineticEnergy(const G4double aValue)
   { fKineticEnergy += aValue; }
     // Kinetic Energy of the track

   inline G4double GetVelocity() const
   { 
     if(fMass==0.){
        return c_light;
     } 
     else{ 
        G4double tMomentum = sqrt(fKineticEnergy*fKineticEnergy +
                                  2.0*fKineticEnergy*fMass);
        G4double tTotalEnergy = fKineticEnergy + fMass; 
        return tMomentum/tTotalEnergy*c_light;
     }   
   }
   // This velocity is the velocity as if in vacuum.
   // (So it is not corrected for the refraction index
   //   in the case of photons - optical or X-rays.)
   // In order to get the velocity in the material, use
   //   GetVelocity of G4Track.
   //

  inline G4double GetBeta() const
   { return (fMass==0.) ? 
      1.0 : 
      sqrt(fKineticEnergy*fKineticEnergy + 2.0*fKineticEnergy*fMass)
      /(fKineticEnergy+fMass); }
    // Velocity of the track in unit of c(light velocity)

   inline G4double GetGamma() const
   { return (fMass==0.) ? DBL_MAX : (fKineticEnergy+fMass)/fMass; }
     // Gamma factor (1/sqrt[1-beta*beta]) of the track    

   inline G4VPhysicalVolume* GetPhysicalVolume()
   { return fpTouchable->GetVolume(); }

   inline G4VTouchable* GetTouchable() const
   { return fpTouchable; }
   inline void SetTouchable(G4VTouchable* apValue)
   { fpTouchable = apValue; }

   inline G4double GetSafety() const
   { return fSafety; }
   inline void SetSafety(const G4double aValue)
   { fSafety = aValue; }

   inline const G4ThreeVector& GetPolarization() const
   { return fPolarization; }
   inline void SetPolarization(const G4ThreeVector& aValue)
   { fPolarization = aValue; }
   inline void AddPolarization(const G4ThreeVector& aValue)
   { fPolarization += aValue; }

   inline G4StepStatus GetStepStatus() const
   { return fStepStatus; }
   inline void SetStepStatus(const G4StepStatus aValue)
   { fStepStatus = aValue; }

   inline const G4VProcess* GetProcessDefinedStep() const
   { return fpProcessDefinedStep; }
     // If the pointer is 0, this means the Step is defined
     // by the user defined limit in the current volume.
   inline void SetProcessDefinedStep(G4VProcess* aValue)
   { fpProcessDefinedStep = aValue; }

   inline G4double GetMass() const
   { return fMass; }
   inline void SetMass(G4double value)
   { fMass = value; }

   inline G4double GetCharge() const
   { return fCharge; }
   inline void SetCharge(G4double value)
   { fCharge = value; }

   inline G4Material* GetMaterial()
   { return fpTouchable->GetVolume()->GetLogicalVolume()->GetMaterial(); }

   inline void SetWeight(G4double aValue)
   { fWeight = aValue; }
   inline G4double GetWeight() const
   { return fWeight; }



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
   G4VTouchable* fpTouchable;
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
};


#endif

