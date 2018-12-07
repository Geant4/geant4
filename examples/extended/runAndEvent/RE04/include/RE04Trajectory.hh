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
/// \file runAndEvent/RE04/include/RE04Trajectory.hh
/// \brief Definition of the RE04Trajectory class
//
//
#ifndef RE04Trajectory_h
#define RE04Trajectory_h 1

#include "G4VTrajectory.hh"
#include "G4Allocator.hh"
#include <stdlib.h>                 // Include from 'system'
#include "G4ios.hh"                 // Include from 'system'
#include <vector>                   // G4RWTValOrderedVector
#include "globals.hh"               // Include from 'global'
#include "G4ParticleDefinition.hh"  // Include from 'particle+matter'
#include "RE04TrajectoryPoint.hh"   // Include from 'tracking'
#include "G4Track.hh"
#include "G4Step.hh"

class G4Polyline;                   // Forward declaration.

typedef std::vector<G4VTrajectoryPoint*>  TrajectoryPointContainer;

//
/// User trajectory class
///
/// - new, delete and "==" operators are overwritten
///
/// - get functions
///     G4int GetTrackID() const, G4int GetParentID() const,
///     G4String GetParticleName() const, G4double GetCharge() const,
///     G4int GetPDGEncoding() const, G4double GetInitialKineticEnergy() const
///     and G4ThreeVector GetInitialMomentum() const
///
/// - void ShowTrajectory(std::ostream& os=G4cout) const
///     invokes the default implementation
///
/// - void DrawTrajectory() const
///     invokes the default implementation
///
/// - void AppendStep(const G4Step* aStep)
///     adds a user trajectory point object, RE04TrajectoryPoint
///
/// - int GetPointEntries() const
///     returns the number of point entries
///
/// - G4VTrajectoryPoint* GetPoint(G4int i) const 
///     gets the i-th trajectory point
///
/// - void MergeTrajectory(G4VTrajectory* secondTrajectory)
///     adds a trajectory to a TrajectoryPointContainer, fPositionRecord
///
/// - G4ParticleDefinition* GetParticleDefinition()
///     get a particle definition from G4ParticleTable
///
/// - const std::map<G4String,G4AttDef>* GetAttDefs() const
///    defines the track ID, the parent ID, the particle name, the charge,
///    the PDG encoding, the initial kinetic energy, the initial momentum,
///    the initial momentum magnitude and the number of points as attiributes 
///
/// - std::vector<G4AttValue>* CreateAttValues() const
///    sets and returns the attributes
//
///////////////////
class RE04Trajectory : public G4VTrajectory
///////////////////
{

//--------
public: 
//--------

// Constructor/Destrcutor

   RE04Trajectory();

   RE04Trajectory(const G4Track* aTrack);
   RE04Trajectory(RE04Trajectory &);
   virtual ~RE04Trajectory();

// Operators
   inline void* operator new(size_t);
   inline void  operator delete(void*);
   inline int operator == (const RE04Trajectory& right) const
   {return (this==&right);} 

// Get/Set functions 
   inline virtual G4int GetTrackID() const
   { return fTrackID; }
   inline virtual G4int GetParentID() const
   { return fParentID; }
   inline virtual G4String GetParticleName() const
   { return fParticleName; }
   inline virtual G4double GetCharge() const
   { return fPDGCharge; }
   inline virtual G4int GetPDGEncoding() const
   { return fPDGEncoding; }
   inline virtual G4double GetInitialKineticEnergy() const
   { return fInitialKineticEnergy; }
   inline virtual G4ThreeVector GetInitialMomentum() const
   { return fInitialMomentum; }

// Other member functions
   virtual void ShowTrajectory(std::ostream& os=G4cout) const;
   //virtual void DrawTrajectory() const;
   virtual void DrawTrajectory() const;
   virtual void AppendStep(const G4Step* aStep);
   virtual int GetPointEntries() const { return fPositionRecord->size(); }
   virtual G4VTrajectoryPoint* GetPoint(G4int i) const 
   { return (*fPositionRecord)[i]; }
   virtual void MergeTrajectory(G4VTrajectory* secondTrajectory);

   G4ParticleDefinition* GetParticleDefinition();

   virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
   virtual std::vector<G4AttValue>* CreateAttValues() const;

//---------
   private:
//---------

  TrajectoryPointContainer* fPositionRecord;
  G4int                     fTrackID;
  G4int                     fParentID;
  G4int                     fPDGEncoding;
  G4double                  fPDGCharge;
  G4String                  fParticleName;
  G4double                  fInitialKineticEnergy;
  G4ThreeVector             fInitialMomentum;

};

extern G4ThreadLocal G4Allocator<RE04Trajectory> * faTrajAllocator;

inline void* RE04Trajectory::operator new(size_t)
{
  if(!faTrajAllocator) faTrajAllocator = new G4Allocator<RE04Trajectory>;
  return (void*)faTrajAllocator->MallocSingle();
}

inline void RE04Trajectory::operator delete(void* aTrajectory)
{
  faTrajAllocator->FreeSingle((RE04Trajectory*)aTrajectory);
}

#endif

