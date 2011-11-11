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


#ifndef RE04Trajectory_h
#define RE04Trajectory_h 1

#include "G4VTrajectory.hh"
#include "G4Allocator.hh"
#include <stdlib.h>                 // Include from 'system'
#include "G4ios.hh"               // Include from 'system'
#include <vector>            // G4RWTValOrderedVector
#include "globals.hh"               // Include from 'global'
#include "G4ParticleDefinition.hh"  // Include from 'particle+matter'
#include "RE04TrajectoryPoint.hh"     // Include from 'tracking'
#include "G4Track.hh"
#include "G4Step.hh"

class G4Polyline;                   // Forward declaration.

typedef std::vector<G4VTrajectoryPoint*>  TrajectoryPointContainer;
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
   inline G4int GetTrackID() const
   { return fTrackID; }
   inline G4int GetParentID() const
   { return fParentID; }
   inline G4String GetParticleName() const
   { return ParticleName; }
   inline G4double GetCharge() const
   { return PDGCharge; }
   inline G4int GetPDGEncoding() const
   { return PDGEncoding; }
   inline G4double GetInitialKineticEnergy() const
   { return initialKineticEnergy; }
   inline G4ThreeVector GetInitialMomentum() const
   { return initialMomentum; }

// Other member functions
   virtual void ShowTrajectory(std::ostream& os=G4cout) const;
   //virtual void DrawTrajectory() const;
   virtual void DrawTrajectory(G4int i_mode = 0) const;
   virtual void AppendStep(const G4Step* aStep);
   virtual int GetPointEntries() const { return positionRecord->size(); }
   virtual G4VTrajectoryPoint* GetPoint(G4int i) const 
   { return (*positionRecord)[i]; }
   virtual void MergeTrajectory(G4VTrajectory* secondTrajectory);

   G4ParticleDefinition* GetParticleDefinition();

   virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
   virtual std::vector<G4AttValue>* CreateAttValues() const;

//---------
   private:
//---------

  TrajectoryPointContainer* positionRecord;
  G4int                     fTrackID;
  G4int                     fParentID;
  G4int                     PDGEncoding;
  G4double                  PDGCharge;
  G4String                  ParticleName;
  G4double                  initialKineticEnergy;
  G4ThreeVector             initialMomentum;

};

extern G4Allocator<RE04Trajectory> aTrajAllocator;

inline void* RE04Trajectory::operator new(size_t)
{
  void* aTrajectory;
  aTrajectory = (void*)aTrajAllocator.MallocSingle();
  return aTrajectory;
}

inline void RE04Trajectory::operator delete(void* aTrajectory)
{
  aTrajAllocator.FreeSingle((RE04Trajectory*)aTrajectory);
}

#endif

