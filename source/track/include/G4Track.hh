// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Track.hh,v 1.1 1999-01-07 16:14:23 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// G4Track.hh
//
// Description:
//   This class represents the partilce under tracking.
//   It includes information related to tracking for examples:
//     1) current position/time of the particle,
//     2) static particle information,
//     3) the pointer to the physical volume where currently
//        the particle exists,
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

#ifndef G4Track_h
#define G4Track_h 1

#include "globals.hh"                 // Include from 'global'
#include "G4ThreeVector.hh"           // Include from 'geometry'
#include "G4LogicalVolume.hh"         // Include from 'geometry'
#include "G4VPhysicalVolume.hh"       // Include from 'geometry'
#include "G4Allocator.hh"             // Include from 'particle+matter'
#include "G4DynamicParticle.hh"       // Include from 'particle+matter'
#include "G4TrackStatus.hh"           // Include from 'tracking'
#include "G4VTouchable.hh"            // Include from 'geometry'

#include "G4Material.hh"

class G4Step;                         // Forward declaration

//////////////
class G4Track
////////////// 
{

//--------
   public:
//--------

// Constructor/Destrcutor
   G4Track();
   G4Track(G4DynamicParticle* apValueDynamicParticle,
           G4double aValueTime,
           const G4ThreeVector& aValuePosition);
      // aValueTime is a global time

   ~G4Track();

// Operators
   inline void *operator new(size_t);
      // Override "new" for "G4Allocator".
   inline void operator delete(void *aTrack);
      // Override "delete" for "G4Allocator".

   int operator==( const G4Track& s);
      // Define "==" operator because "G4TrackVector" uses 
      //"RWPtrOrderdVector" which requires this.

// Get/Set functions
   const G4ThreeVector& GetPosition() const;
   void SetPosition(const G4ThreeVector& aValue);

   G4double GetGlobalTime() const;
   void SetGlobalTime(const G4double aValue);
     // Time since the event in which the track belongs is created.

   G4double GetLocalTime() const;
   void SetLocalTime(const G4double aValue);
      // Time since the current track is created.

   G4double GetProperTime() const;
   void SetProperTime(const G4double aValue);
      // Proper time of the current track

   G4double GetTrackLength() const;
   void AddTrackLength(const G4double aValue);
      // Accumulated the track length

   G4int GetParentID() const;
   void SetParentID(const G4int aValue);

   G4int GetTrackID() const;
   void SetTrackID(const G4int aValue);

   G4VPhysicalVolume* GetVolume() const;

   G4VPhysicalVolume* GetNextVolume() const;

   G4Material* GetMaterial() const;
   G4Material* GetNextMaterial() const;

   G4VTouchable* GetTouchable() const;
   void SetTouchable(G4VTouchable* apValue);

   G4VTouchable* GetNextTouchable() const;
   void SetNextTouchable(G4VTouchable* apValue);

   G4double GetKineticEnergy() const;
   void SetKineticEnergy(const G4double aValue);

   G4double GetVelocity() const;

   const G4ThreeVector& GetMomentumDirection() const;
   void SetMomentumDirection(const G4ThreeVector& aValue);

   const G4ThreeVector& GetPolarization() const;
   void SetPolarization(const G4ThreeVector& aValue);

   G4TrackStatus GetTrackStatus() const;
   void SetTrackStatus(const G4TrackStatus aTrackStatus);

   G4bool IsBelowThreshold() const;
   void   SetBelowThresholdFlag(G4bool value = true);

   G4bool IsGoodForTracking() const;
   void   SetGoodForTrackingFlag(G4bool value = true);

   G4int GetCurrentStepNumber() const;
   void IncrementCurrentStepNumber();

   G4double GetTotalEnergy() const;

   G4ThreeVector GetMomentum() const;

   const G4DynamicParticle* GetDynamicParticle() const;

   G4ParticleDefinition* GetDefinition() const;

   G4double GetStepLength() const;
   void SetStepLength(G4double value);

   G4Step* GetStep() const;
   void SetStep(G4Step* aValue);

   const G4ThreeVector& GetVertexPosition() const;
   void SetVertexPosition(const G4ThreeVector& aValue);

   const G4ThreeVector& GetVertexMomentumDirection() const;
   void SetVertexMomentumDirection(const G4ThreeVector& aValue);

   G4double GetVertexKineticEnergy() const;
   void SetVertexKineticEnergy(const G4double aValue);

   G4LogicalVolume* GetLogicalVolumeAtVertex() const;
   void SetLogicalVolumeAtVertex(G4LogicalVolume* );

   const G4VProcess* GetCreatorProcess() const;
   void SetCreatorProcess(G4VProcess* aValue);

   G4double GetWeight() const;
   void     SetWeight(G4double aValue);

//---------
   private:
//---------

// Member data
   G4int fCurrentStepNumber;       // Total steps number up to now
   G4ThreeVector fPosition;        // Current positon
   G4double fGlobalTime;           // Time since the event is created
   G4double fLocalTime;            // Time since the track is created
   G4double fTrackLength;          // Accumulated track length
   G4int fParentID;
   G4int fTrackID;

   G4VTouchable* fpTouchable;
   G4VTouchable* fpNextTouchable;

   G4DynamicParticle* fpDynamicParticle;
   G4TrackStatus fTrackStatus;

   G4bool  fBelowThreshold;
   // This flag is set to true if this track energy is below
   // threshold energy in this material determined by the range cut value
   G4bool  fGoodForTracking;
   // This flag is set by processes if this track should be tracked
   // even if the energy is below threshold

   G4double fStepLength;           
      // Before the end of the AlongStepDoIt loop, this keeps the initial 
      // Step length which is determined by the shortest geometrical Step 
      // proposed by a physics process. After finishing the AlongStepDoIt,
      // this will be set equal to 'StepLength' in G4Step.

   G4double fWeight;
     // This is a weight for this track used by G4VEvtBiasMechanism 
     // to execute inclusive simulation for hadronic/electomagnetic shower
     // and neutron transportation etc. 

   G4Step* fpStep;

   G4ThreeVector fVtxPosition;          // (x,y,z) of the vertex
   G4ThreeVector fVtxMomentumDirection; // Momentum direction at the vertex
   G4double fVtxKineticEnergy;          // Kinetic energy at the vertex
   G4LogicalVolume* fpLVAtVertex;       //Logical Volume at the vertex
   G4VProcess* fpCreatorProcess;        // Process which created the track
   
};


//-----------------------------------------------------------------
// Definitions of inline functions
//-----------------------------------------------------------------

// Operators
   extern G4Allocator<G4Track> aTrackAllocator;
   inline void* G4Track::operator new(size_t)
   { void *aTrack;
     aTrack = (void *) aTrackAllocator.MallocSingle();
     return aTrack;
   }
      // Override "new" for "G4Allocator".

   inline void G4Track::operator delete(void *aTrack)
   { aTrackAllocator.FreeSingle((G4Track *) aTrack);}
      // Override "delete" for "G4Allocator".

   inline int G4Track::operator==( const G4Track& s)
   { return (this==&s) ? 1 : 0; }
      // Define "==" operator because "G4TrackVector" uses 
      // "RWPtrOrderdVector" which requires this.

// Get/Set functions
   inline const G4ThreeVector& G4Track::GetPosition() const
   { return fPosition; }
   inline void G4Track::SetPosition(const G4ThreeVector& aValue)
   { fPosition = aValue; }

   inline G4double G4Track::GetGlobalTime() const
   { return fGlobalTime; }
   inline void G4Track::SetGlobalTime(const G4double aValue)
   { fGlobalTime = aValue; }
     // Time since the event in which the track belongs is created.

   inline G4double G4Track::GetLocalTime() const
   { return fLocalTime; }
   inline void G4Track::SetLocalTime(const G4double aValue)
   { fLocalTime = aValue; }
      // Time since the current track is created.

   inline G4double G4Track::GetProperTime() const
   { return fpDynamicParticle->GetProperTime(); }
   inline void G4Track::SetProperTime(const G4double aValue)
   { fpDynamicParticle->SetProperTime(aValue); }
      // Proper time of the current track

   inline G4double G4Track::GetTrackLength() const
   { return fTrackLength; }
   inline void G4Track::AddTrackLength(const G4double aValue)
   { fTrackLength += aValue; }
      // Accumulated track length

   inline G4int G4Track::GetParentID() const
   { return fParentID; }
   inline void G4Track::SetParentID(const G4int aValue)
   { fParentID = aValue; }

   inline G4int G4Track::GetTrackID() const
   { return fTrackID; }
   inline void G4Track::SetTrackID(const G4int aValue)
   { fTrackID = aValue; }

   inline G4VPhysicalVolume* G4Track::GetVolume() const
   { return fpTouchable->GetVolume(); }

   inline G4VPhysicalVolume* G4Track::GetNextVolume() const
   { return fpNextTouchable->GetVolume(); }

   inline G4Material* G4Track::GetMaterial() const
   { return fpTouchable->GetVolume()->GetLogicalVolume()->GetMaterial(); }

   inline G4Material* G4Track::GetNextMaterial() const
   { return fpNextTouchable->GetVolume()->GetLogicalVolume()->GetMaterial(); }

   inline G4VTouchable* G4Track::GetTouchable() const
   { return fpTouchable; }
   inline void G4Track::SetTouchable(G4VTouchable* apValue)
   { fpTouchable = apValue; }

   inline G4VTouchable* G4Track::GetNextTouchable() const
   { return fpNextTouchable; }
   inline void G4Track::SetNextTouchable(G4VTouchable* apValue)
   { fpNextTouchable = apValue; }

   inline G4double G4Track::GetKineticEnergy() const
   { return fpDynamicParticle->GetKineticEnergy(); }
   inline void G4Track::SetKineticEnergy(const G4double aValue)
   { fpDynamicParticle->SetKineticEnergy(aValue); }

   inline G4double G4Track::GetVelocity() const
   { 
    G4double velocity ;

    G4double mass = fpDynamicParticle->GetMass();
    if( mass == 0. )
    {
     velocity = c_light ; 
     if((fpDynamicParticle->GetDefinition()->GetParticleName() ==
                                                            "gamma")
        ||
        (fpDynamicParticle->GetDefinition()->GetParticleName() ==
                                                    "opticalphoton"))
     {
       G4Material*
        mat=fpTouchable->GetVolume()->GetLogicalVolume()->GetMaterial();
       if(mat->GetMaterialPropertiesTable() != NULL)
       {
        if(mat->GetMaterialPropertiesTable()->GetProperty("RINDEX") != NULL ) 
          velocity /= 
          mat->GetMaterialPropertiesTable()->GetProperty("RINDEX")->
          GetMinProperty() ; 
       }
     }  
       
    }
    else
    {
     G4double T = fpDynamicParticle->GetKineticEnergy();
     velocity = c_light*sqrt(T*(T+2.*mass))/(T+mass) ;
    }

    return velocity ;

   }

   inline const G4ThreeVector& G4Track::GetMomentumDirection() const
   { return fpDynamicParticle->GetMomentumDirection(); }
   inline void G4Track::SetMomentumDirection(const G4ThreeVector& aValue)
   { fpDynamicParticle->SetMomentumDirection(aValue) ;}

   inline const G4ThreeVector& G4Track::GetPolarization() const
   { return fpDynamicParticle->GetPolarization(); }
   inline void G4Track::SetPolarization(const G4ThreeVector& aValue)
   { fpDynamicParticle->SetPolarization(aValue.x(),
                                        aValue.y(),
                                        aValue.z()); }

   inline G4TrackStatus G4Track::GetTrackStatus() const
   { return fTrackStatus; }
   inline void G4Track::SetTrackStatus(const G4TrackStatus aTrackStatus)
   { fTrackStatus = aTrackStatus; }

   inline G4int G4Track::GetCurrentStepNumber() const
   { return fCurrentStepNumber; }
   inline void G4Track::IncrementCurrentStepNumber()
   { fCurrentStepNumber++; }

   inline G4double G4Track::GetTotalEnergy() const
   { return fpDynamicParticle->GetTotalEnergy(); }

   inline G4ThreeVector G4Track::GetMomentum() const
   { return fpDynamicParticle->GetMomentum(); }

   inline const G4DynamicParticle* G4Track::GetDynamicParticle() const
   { return fpDynamicParticle; }

   inline G4ParticleDefinition* G4Track::GetDefinition() const
   { return fpDynamicParticle->GetDefinition(); }

   inline G4double G4Track::GetStepLength() const
   { return fStepLength; }
   inline void G4Track::SetStepLength(G4double value)
   { fStepLength = value; }

   inline const G4ThreeVector& G4Track::GetVertexPosition() const
   { return fVtxPosition; }
   inline void G4Track::SetVertexPosition(const G4ThreeVector& aValue)
   { fVtxPosition = aValue; }

   inline const G4ThreeVector& G4Track::GetVertexMomentumDirection() const
   { return fVtxMomentumDirection; }
   inline void G4Track::SetVertexMomentumDirection(const G4ThreeVector& aValue)
   { fVtxMomentumDirection = aValue ;}

   inline G4double G4Track::GetVertexKineticEnergy() const
   { return fVtxKineticEnergy; }
   inline void G4Track::SetVertexKineticEnergy(const G4double aValue)
   { fVtxKineticEnergy = aValue; }

   inline  G4LogicalVolume* G4Track::GetLogicalVolumeAtVertex() const
   { return fpLVAtVertex; } 
   inline void G4Track::SetLogicalVolumeAtVertex(G4LogicalVolume* aValue)
   { fpLVAtVertex = aValue; }

   inline const G4VProcess* G4Track::GetCreatorProcess() const
   { return fpCreatorProcess; }
     // If the pointer is NULL, this means the track is created
     // by the event generator, i.e. the primary track.If it is not
     // NULL, it points to the process which created this track.
   inline void G4Track::SetCreatorProcess(G4VProcess* aValue)
   { fpCreatorProcess = aValue; }

   inline G4bool G4Track::IsBelowThreshold() const
   { return fBelowThreshold; }
   inline void    G4Track::SetBelowThresholdFlag(G4bool value)
   { fBelowThreshold = value; }

   inline G4bool  G4Track::IsGoodForTracking() const
   { return fGoodForTracking; }
   inline void    G4Track::SetGoodForTrackingFlag(G4bool value)
   { fGoodForTracking = value; }

   inline void  G4Track::SetWeight(G4double aValue)
   { fWeight = aValue; }
   inline G4double G4Track::GetWeight() const
   { return fWeight; }





//-------------------------------------------------------------
// To implement bi-directional association between G4Step and
// and G4Track, a combined usage of 'forward declaration' and
// 'include' is necessary.
//-------------------------------------------------------------
#include "G4Step.hh"

   inline G4Step* G4Track::GetStep() const
   { return fpStep; }
   inline void G4Track::SetStep(G4Step* aValue)
   { fpStep = aValue; }

#endif









