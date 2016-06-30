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
// $Id: G4Track.hh 94983 2016-01-13 11:02:33Z gcosmo $
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
//        the particle exists
//
//---------------------------------------------------------------
//   Modification for G4TouchableHandle             22 Oct. 2001  R.Chytracek
//   Add MaterialCutCouple                          08 Oct. 2002  H.Kurashige
//   Add SetVelocityTableProperties                 02 Apr. 2011  H.Kurashige
//   Add fVelocity and Set/GetVelocity              29 Apr. 2011  H.Kurashige
//   Use G4VelocityTable                     17 AUg. 2011 H.Kurashige

#ifndef G4Track_h
#define G4Track_h 1

#include <cmath>                      // Include from 'system'

#include "globals.hh"                 // Include from 'global'
#include "trkdefs.hh"                 // Include DLL defs...
#include "G4ThreeVector.hh"           // Include from 'geometry'
#include "G4LogicalVolume.hh"         // Include from 'geometry'
#include "G4VPhysicalVolume.hh"       // Include from 'geometry'
#include "G4Allocator.hh"             // Include from 'particle+matter'
#include "G4DynamicParticle.hh"       // Include from 'particle+matter'
#include "G4TrackStatus.hh"           // Include from 'tracking'
#include "G4TouchableHandle.hh"       // Include from 'geometry'
#include "G4VUserTrackInformation.hh"
#include "G4PhysicsModelCatalog.hh"

#include "G4Material.hh"

class G4Step;                         // Forward declaration
class G4MaterialCutsCouple;
class G4VelocityTable;
class G4VAuxiliaryTrackInformation;

#include <map>

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
   G4Track(const G4Track&);
   // Copy Constructor copys members other than tracking information
 
private:
  // Hide assignment operator as private
   G4Track& operator=(const  G4Track&);

//--------
public: // With description

// Destrcutor
   ~G4Track();

// Operators
   inline void *operator new(size_t);
      // Override "new" for "G4Allocator".
   inline void operator delete(void *aTrack);
      // Override "delete" for "G4Allocator".

   G4bool operator==( const G4Track& );
  
//--------
public: // With description
// Copy information of the track (w/o tracking information)
   void CopyTrackInfo(const G4Track&);

// Get/Set functions
  // track ID
   G4int GetTrackID() const;
   void SetTrackID(const G4int aValue);

   G4int GetParentID() const;
   void SetParentID(const G4int aValue);

  // dynamic particle 
   const G4DynamicParticle* GetDynamicParticle() const;

  // particle definition
    const G4ParticleDefinition* GetParticleDefinition() const;
   //  following method of GetDefinition remains 
   //  because of backward compatiblity. It will be removed in future 
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

   const G4MaterialCutsCouple* GetMaterialCutsCouple() const;
   const G4MaterialCutsCouple* GetNextMaterialCutsCouple() const;

   const G4VTouchable*      GetTouchable() const;
   const G4TouchableHandle& GetTouchableHandle() const;
   void SetTouchableHandle( const G4TouchableHandle& apValue);

   const G4VTouchable*      GetNextTouchable() const;
   const G4TouchableHandle& GetNextTouchableHandle() const;
   void SetNextTouchableHandle( const G4TouchableHandle& apValue);

   const G4VTouchable*      GetOriginTouchable() const;
   const G4TouchableHandle& GetOriginTouchableHandle() const;
   void SetOriginTouchableHandle( const G4TouchableHandle& apValue);

  // energy
   G4double GetKineticEnergy() const;
   void SetKineticEnergy(const G4double aValue);

   G4double GetTotalEnergy() const;

 
  // moemtnum
   const G4ThreeVector& GetMomentumDirection() const;
   void SetMomentumDirection(const G4ThreeVector& aValue);

   G4ThreeVector GetMomentum() const;

   // velocity
   G4double GetVelocity() const;
   void     SetVelocity(G4double val);
 
   G4double CalculateVelocity() const;
   G4double CalculateVelocityForOpticalPhoton() const;

   G4bool   UseGivenVelocity() const; 
   void     UseGivenVelocity(G4bool val);

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

   const G4LogicalVolume* GetLogicalVolumeAtVertex() const;
   void SetLogicalVolumeAtVertex(const G4LogicalVolume* );

   const G4VProcess* GetCreatorProcess() const;
   void SetCreatorProcess(const G4VProcess* aValue);

   inline void SetCreatorModelIndex(G4int idx);

   inline const G4String& GetCreatorModelName() const;

   inline G4int GetCreatorModelID() const;

  // track weight
  // These are methods for manipulating a weight for this track.
   G4double GetWeight() const;
   void     SetWeight(G4double aValue);

  // User information
  G4VUserTrackInformation* GetUserInformation() const;
  void SetUserInformation(G4VUserTrackInformation* aValue) const;
 
  // Velocity table
  static void SetVelocityTableProperties(G4double t_max, G4double t_min, G4int nbin);
  static G4double GetMaxTOfVelocityTable();
  static G4double GetMinTOfVelocityTable();
  static G4int    GetNbinOfVelocityTable();

//---------
   private:
//---------
   // Member data
   G4int fCurrentStepNumber;       // Total steps number up to now
   G4ThreeVector fPosition;        // Current positon
   G4double fGlobalTime;           // Time since the event is created
   G4double fLocalTime;            // Time since the track is created
   G4double fTrackLength;          // Accumulated track length
   G4int    fParentID;
   G4int    fTrackID;
   G4double fVelocity; 

   G4TouchableHandle fpTouchable;
   G4TouchableHandle fpNextTouchable;
   G4TouchableHandle fpOriginTouchable;
  // Touchable Handle

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
     // This is a weight for this track 

   const G4Step* fpStep;

   G4ThreeVector fVtxPosition;          // (x,y,z) of the vertex
   G4ThreeVector fVtxMomentumDirection; // Momentum direction at the vertex
   G4double fVtxKineticEnergy;          // Kinetic energy at the vertex
   const G4LogicalVolume* fpLVAtVertex; //Logical Volume at the vertex
   const G4VProcess* fpCreatorProcess; // Process which created the track
   G4int fCreatorModelIndex;           // Index of the physics model which created the track
   
   mutable G4VUserTrackInformation* fpUserInformation;

   // cached values for CalculateVelocity  
   mutable G4Material*               prev_mat;
   mutable G4MaterialPropertyVector* groupvel;
   mutable G4double                  prev_velocity;
   mutable G4double                  prev_momentum;

   G4bool          is_OpticalPhoton; 
   static G4ThreadLocal G4VelocityTable*  velTable;
 
   G4bool          useGivenVelocity;
      // do not calclulate velocity and just use current fVelocity
      // if this flag is set

   mutable std::map<G4int,G4VAuxiliaryTrackInformation*>* fpAuxiliaryTrackInformationMap;

//--------
public: 
//--------

   void SetAuxiliaryTrackInformation(G4int idx, G4VAuxiliaryTrackInformation* info) const;
   G4VAuxiliaryTrackInformation* GetAuxiliaryTrackInformation(G4int idx) const;
   std::map<G4int,G4VAuxiliaryTrackInformation*>* GetAuxiliaryTrackInformationMap() const
   { return fpAuxiliaryTrackInformationMap; }

   void RemoveAuxiliaryTrackInformation(G4int idx);
   void RemoveAuxiliaryTrackInformation(G4String& name);
      // Note: G4VAuxiliaryTrackInformation object itself is *NOT* deleted

//--------
private:
//--------

  void ClearAuxiliaryTrackInformation();
};

#include "G4Track.icc"

#endif
