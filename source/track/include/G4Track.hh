// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Track.hh,v 1.9 2000-10-18 15:00:09 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// G4Track.hh
//
// Class Description:
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
#include "G4VUserTrackInformation.hh"

#include "G4Material.hh"

class G4Step;                         // Forward declaration

//////////////
class G4Track
////////////// 
{

//--------
public: // With description

// Constructor
   G4Track();
   G4Track(G4DynamicParticle* apValueDynamicParticle,
           G4double aValueTime,
           const G4ThreeVector& aValuePosition);
      // aValueTime is a global time

//--------
public: 

// Destrcutor
   ~G4Track();

// Operators
   inline void *operator new(size_t);
      // Override "new" for "G4Allocator".
   inline void operator delete(void *aTrack);
      // Override "delete" for "G4Allocator".

   int operator==( const G4Track& s);
      // Define "==" operator because "G4TrackVector" uses 
      //"RWPtrOrderdVector" which requires this.

//--------
public: // With description

// Get/Set functions
  // track ID
   G4int GetTrackID() const;
   void SetTrackID(const G4int aValue);

   G4int GetParentID() const;
   void SetParentID(const G4int aValue);

  // dynamic particle 
   const G4DynamicParticle* GetDynamicParticle() const;

  // particle definition
   G4ParticleDefinition* GetDefinition() const;

   // position, time 
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

  // volume, material, touchable
   G4VPhysicalVolume* GetVolume() const;
   G4VPhysicalVolume* GetNextVolume() const;

   G4Material* GetMaterial() const;
   G4Material* GetNextMaterial() const;

   const G4VTouchable* GetTouchable() const;
   void SetTouchable(const G4VTouchable* apValue);

   const G4VTouchable* GetNextTouchable() const;
   void SetNextTouchable(const G4VTouchable* apValue);

  // energy
   G4double GetKineticEnergy() const;
   void SetKineticEnergy(const G4double aValue);

   G4double GetTotalEnergy() const;

 
  // moemtnum
   const G4ThreeVector& GetMomentumDirection() const;
   void SetMomentumDirection(const G4ThreeVector& aValue);

   G4ThreeVector GetMomentum() const;

   G4double GetVelocity() const;


  // polarization 
   const G4ThreeVector& GetPolarization() const;
   void SetPolarization(const G4ThreeVector& aValue);

  // track status, flags for tracking
   G4TrackStatus GetTrackStatus() const;
   void SetTrackStatus(const G4TrackStatus aTrackStatus);

   G4bool IsBelowThreshold() const;
   void   SetBelowThresholdFlag(G4bool value = true);
     // The flag of "BelowThreshold" is set to true 
     // if this track energy is below threshold energy 
     //  in this material determined by the range cut value

   G4bool IsGoodForTracking() const;
   void   SetGoodForTrackingFlag(G4bool value = true);
     // The flag of "GoodForTracking" is set by processes 
     // if this track should be tracked
     // even if the energy is below threshold

  // track length
   G4double GetTrackLength() const;
   void AddTrackLength(const G4double aValue);
      // Accumulated the track length

  // step information
   const G4Step* GetStep() const;
   void SetStep(const G4Step* aValue);

   G4int GetCurrentStepNumber() const;
   void IncrementCurrentStepNumber();

   G4double GetStepLength() const;
   void SetStepLength(G4double value);
      // Before the end of the AlongStepDoIt loop,StepLength keeps 
      // the initial value which is determined by the shortest geometrical Step 
      // proposed by a physics process. After finishing the AlongStepDoIt,
      // it will be set equal to 'StepLength' in G4Step.

  // vertex (,where this track was created) information  
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

  // track weight
     // These are methods for manipulating a weight for this track.
     // The track weight is used by G4VEvtBiasMechanism 
     // to execute inclusive simulation for hadronic/electomagnetic shower
     // and neutron transportation etc. 
   G4double GetWeight() const;
   void     SetWeight(G4double aValue);

  // User information
  G4VUserTrackInformation* GetUserInformation() const;
  void SetUserInformation(G4VUserTrackInformation* aValue);

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

   const G4VTouchable* fpTouchable;
   const G4VTouchable* fpNextTouchable;

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

   const G4Step* fpStep;

   G4ThreeVector fVtxPosition;          // (x,y,z) of the vertex
   G4ThreeVector fVtxMomentumDirection; // Momentum direction at the vertex
   G4double fVtxKineticEnergy;          // Kinetic energy at the vertex
   G4LogicalVolume* fpLVAtVertex;       //Logical Volume at the vertex
   G4VProcess* fpCreatorProcess;        // Process which created the track
   
   G4VUserTrackInformation* fpUserInformation;
};
#include "G4Step.hh"
#include "G4Track.icc"

#endif















