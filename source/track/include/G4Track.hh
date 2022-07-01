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
// G4Track
//
// Class description:
//
// This class describes the particle under tracking.
// It includes information related to tracking, i.e.:
//     1) current position/time of the particle,
//     2) static particle information,
//     3) the pointer to the physical volume where currently
//        the particle exists

// Author: Katsuya Amako, KEK - 1995
// Revisions: Hisaya Kurashige, 1998-2011
// --------------------------------------------------------------------
#ifndef G4Track_hh
#define G4Track_hh 1

#include <cmath>  // Include from 'system'
#include <map>
#include <CLHEP/Units/PhysicalConstants.h>

#include "globals.hh"            // Include from 'global'
#include "trkdefs.hh"            // Include DLL defs...
#include "G4ThreeVector.hh"      // Include from 'geometry'
#include "G4LogicalVolume.hh"    // Include from 'geometry'
#include "G4VPhysicalVolume.hh"  // Include from 'geometry'
#include "G4Allocator.hh"        // Include from 'particle+matter'
#include "G4DynamicParticle.hh"  // Include from 'particle+matter'
#include "G4TrackStatus.hh"      // Include from 'tracking'
#include "G4TouchableHandle.hh"  // Include from 'geometry'
#include "G4VUserTrackInformation.hh"
#include "G4PhysicsModelCatalog.hh"
#include "G4Material.hh"
#include "G4Profiler.hh"

class G4Step;  // Forward declaration
class G4MaterialCutsCouple;
class G4VAuxiliaryTrackInformation;
class G4VProcess;

class G4Track
{
  public:
  // the profiler aliases are only used when compiled with GEANT4_USE_TIMEMORY
  using ProfilerConfig = G4ProfilerConfig<G4ProfileType::Track>;

    G4Track();
      // Default constructor
    G4Track(G4DynamicParticle* apValueDynamicParticle,
            G4double aValueTime,
            const G4ThreeVector& aValuePosition);
      // Constructor - aValueTime is a global time

    G4Track(const G4Track&);
      // Copy Constructor - copies members other than tracking information

    ~G4Track();
      // Destructor

    inline void* operator new(std::size_t);
      // Override "new" for "G4Allocator".
    inline void operator delete(void* aTrack);
      // Override "delete" for "G4Allocator".

    G4Track& operator=(const G4Track&);
      // Assignment operator

    inline G4bool operator==(const G4Track&);
    inline G4bool operator!=(const G4Track&);
      // Equality operators

    void CopyTrackInfo(const G4Track&);
      // Copy information of the track (w/o tracking information)

    G4int GetTrackID() const;
    void SetTrackID(const G4int aValue);
      // Get/Set functions track ID

    G4int GetParentID() const;
    void SetParentID(const G4int aValue);

    const G4DynamicParticle* GetDynamicParticle() const;
      // Dynamic particle

    const G4ParticleDefinition* GetParticleDefinition() const;
      // Particle definition
    G4ParticleDefinition* GetDefinition() const;
      // Obsolete, for backwards compatibility...

    const G4ThreeVector& GetPosition() const;
    void SetPosition(const G4ThreeVector& aValue);
      // Position, time

    G4double GetGlobalTime() const;
    void SetGlobalTime(const G4double aValue);
      // Time since the event in which the track belongs is created

    G4double GetLocalTime() const;
    void SetLocalTime(const G4double aValue);
      // Time since the current track is created

    G4double GetProperTime() const;
    void SetProperTime(const G4double aValue);
      // Proper time of the current track

    G4VPhysicalVolume* GetVolume() const;
    G4VPhysicalVolume* GetNextVolume() const;
      // Volume, material, touchable

    G4Material* GetMaterial() const;
    G4Material* GetNextMaterial() const;

    const G4MaterialCutsCouple* GetMaterialCutsCouple() const;
    const G4MaterialCutsCouple* GetNextMaterialCutsCouple() const;

    const G4VTouchable* GetTouchable() const;
    const G4TouchableHandle& GetTouchableHandle() const;
    void SetTouchableHandle(const G4TouchableHandle& apValue);

    const G4VTouchable* GetNextTouchable() const;
    const G4TouchableHandle& GetNextTouchableHandle() const;
    void SetNextTouchableHandle(const G4TouchableHandle& apValue);

    const G4VTouchable* GetOriginTouchable() const;
    const G4TouchableHandle& GetOriginTouchableHandle() const;
    void SetOriginTouchableHandle(const G4TouchableHandle& apValue);

    G4double GetKineticEnergy() const;
    void SetKineticEnergy(const G4double aValue);
      // Energy

    G4double GetTotalEnergy() const;

    const G4ThreeVector& GetMomentumDirection() const;
    void SetMomentumDirection(const G4ThreeVector& aValue);
      // Momentum

    G4ThreeVector GetMomentum() const;

    G4double GetVelocity() const;
    void SetVelocity(G4double val);
      // Velocity

    G4double CalculateVelocity() const;
    G4double CalculateVelocityForOpticalPhoton() const;

    G4bool UseGivenVelocity() const;
    void UseGivenVelocity(G4bool val);

    const G4ThreeVector& GetPolarization() const;
    void SetPolarization(const G4ThreeVector& aValue);
      // Polarization

    G4TrackStatus GetTrackStatus() const;
    void SetTrackStatus(const G4TrackStatus aTrackStatus);
      // Track status, flags for tracking

    G4bool IsBelowThreshold() const;
    void SetBelowThresholdFlag(G4bool value = true);
      // The flag of "BelowThreshold" is set to true
      // If this track energy is below threshold energy
      // in this material is determined by the range cut value

    G4bool IsGoodForTracking() const;
    void SetGoodForTrackingFlag(G4bool value = true);
      // The flag of "GoodForTracking" is set by processes
      // if this track should be tracked
      // even if the energy is below threshold

    G4double GetTrackLength() const;
    void AddTrackLength(const G4double aValue);
      // Accumulated track length

    const G4Step* GetStep() const;
    void SetStep(const G4Step* aValue);
      // Step information

    G4int GetCurrentStepNumber() const;
    void IncrementCurrentStepNumber();

    G4double GetStepLength() const;
    void SetStepLength(G4double value);
      // Before the end of the AlongStepDoIt() loop, StepLength keeps
      // the initial value which is determined by the shortest geometrical Step
      // proposed by a physics process. After finishing the AlongStepDoIt(),
      // it will be set equal to 'StepLength' in G4Step

    const G4ThreeVector& GetVertexPosition() const;
    void SetVertexPosition(const G4ThreeVector& aValue);
      // Vertex (where this track was created) information

    const G4ThreeVector& GetVertexMomentumDirection() const;
    void SetVertexMomentumDirection(const G4ThreeVector& aValue);

    G4double GetVertexKineticEnergy() const;
    void SetVertexKineticEnergy(const G4double aValue);

    const G4LogicalVolume* GetLogicalVolumeAtVertex() const;
    void SetLogicalVolumeAtVertex(const G4LogicalVolume*);

    const G4VProcess* GetCreatorProcess() const;
    void SetCreatorProcess(const G4VProcess* aValue);

    inline void SetCreatorModelID(const G4int id);
    inline G4int GetCreatorModelID() const;  
    inline G4int GetCreatorModelIndex() const;
    inline const G4String GetCreatorModelName() const;
      // Identification of the physics model that created the track:
      // each of the three information (ID, index, name) is unique
      // (the model ID and its name are supposed to be used in Geant4
      // code, whereas the index is meant for plotting in user code)

    inline const G4ParticleDefinition* GetParentResonanceDef() const;
    inline void SetParentResonanceDef(const G4ParticleDefinition* parent);
    inline G4int GetParentResonanceID() const;
    inline void SetParentResonanceID(const G4int parentID );
    inline G4bool HasParentResonance() const;
    inline G4int GetParentResonancePDGEncoding() const;
    inline G4String GetParentResonanceName() const;
    inline G4double GetParentResonanceMass() const;
      // Because short-lived resonances (e.g. omega, phi, rho, delta, etc.)
      // do not have corresponding track objects, if the track is produced
      // by a resonance parent, these methods allow to get/set information
      // regarding this short-lived parent.
      // The ID is a unique (integer) identifier for each resonance (which
      // corresponds to the rounded integer of the mass of the resonance
      // in keV), which allows to know if two (or more) tracks originated
      // from the same parent resonance: this should not be confused with
      // the parent-track-ID (fParentID) which corresponds to its closest
      // ancestor which is not a short-lived resonance (and therefore has
      // a corresponding track object).
      // In the case of a track non originating from a resonance parent,
      // the above "Get" methods return, respectively: nullptr, 0, false,
      // 0, "", 0.
  
    G4double GetWeight() const;
    void SetWeight(G4double aValue);
      // Track weight; methods for manipulating a weight for this track

    G4VUserTrackInformation* GetUserInformation() const;
    void SetUserInformation(G4VUserTrackInformation* aValue) const;
      // User information

    void SetAuxiliaryTrackInformation(G4int id,
                                      G4VAuxiliaryTrackInformation* info) const;
    G4VAuxiliaryTrackInformation* GetAuxiliaryTrackInformation(G4int id) const;
    inline std::map<G4int, G4VAuxiliaryTrackInformation*>*
           GetAuxiliaryTrackInformationMap() const;

    void RemoveAuxiliaryTrackInformation(G4int id);
    void RemoveAuxiliaryTrackInformation(G4String& name);
      // Note: G4VAuxiliaryTrackInformation object itself is *NOT* deleted

  private:

    void ClearAuxiliaryTrackInformation();

    // Member data

    G4ThreeVector fPosition;
      // Current positon
    G4double fGlobalTime = 0.0;
      // Time since the event is created
    G4double fLocalTime = 0.0;
      // Time since the track is created
    G4double fTrackLength = 0.0;
      // Accumulated track length
 
    G4double fVelocity = 0.0;

    G4TouchableHandle fpTouchable;
    G4TouchableHandle fpNextTouchable;
    G4TouchableHandle fpOriginTouchable;
      // Touchable Handle

    G4DynamicParticle* fpDynamicParticle = nullptr;
    G4TrackStatus fTrackStatus = fAlive;

    G4double fStepLength = 0.0;
      // Before the end of the AlongStepDoIt loop, this keeps the initial
      // Step length which is determined by the shortest geometrical Step
      // proposed by a physics process. After finishing the AlongStepDoIt,
      // this will be set equal to 'StepLength' in G4Step.

    G4double fWeight = 1.0;
      // This is a weight for this track

    const G4Step* fpStep = nullptr;

    G4ThreeVector fVtxPosition;
      // (x,y,z) of the vertex
    G4ThreeVector fVtxMomentumDirection;
      // Momentum direction at the vertex
    G4double fVtxKineticEnergy = 0.0;
      // Kinetic energy at the vertex
    const G4LogicalVolume* fpLVAtVertex = nullptr;
      // Logical Volume at the vertex
    const G4VProcess* fpCreatorProcess = nullptr;
      // Process which created the track

    mutable G4VUserTrackInformation* fpUserInformation = nullptr;

    mutable G4Material* prev_mat = nullptr;
    mutable G4MaterialPropertyVector* groupvel = nullptr;
    mutable G4double prev_velocity = 0.0;
    mutable G4double prev_momentum = 0.0;
      // cached values for CalculateVelocity

    mutable std::map<G4int, G4VAuxiliaryTrackInformation*>*
            fpAuxiliaryTrackInformationMap = nullptr;

    G4int fCurrentStepNumber = 0;
      // Total steps number up to now

    G4int fCreatorModelID = -1;
      // ID of the physics model which created the track

    const G4ParticleDefinition* fParentResonanceDef = nullptr;
      // Pointer to the particle definition of a short-lived resonance,
      // in the case that the track is produced by a resonance parent
      // (which does not have a corresponding track object)
    G4int fParentResonanceID = 0;
      // Unique ID for the parent resonance, in the case that the track
      // is produced by a resonance parent, else 0

    G4int fParentID = 0;
    G4int fTrackID = 0;

    G4bool fBelowThreshold = false;
      // This flag is set to true if this track energy is below
      // threshold energy in this material determined by the range cut value
    G4bool fGoodForTracking = false;
      // This flag is set by processes if this track should be tracked
      // even if the energy is below threshold

    G4bool is_OpticalPhoton = false;

    G4bool useGivenVelocity = false;
      // do not calculate velocity and just use current fVelocity
      // if this flag is set
};

#include "G4Track.icc"

#endif
