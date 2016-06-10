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
// $Id: G4SmoothTrajectory.hh 69003 2013-04-15 09:25:23Z gcosmo $
//
//---------------------------------------------------------------
//
// G4SmoothTrajectory.hh
//
// class description:
//   This class represents the trajectory of a particle tracked.
//   It includes information of 
//     1) List of trajectory points which compose the trajectory,
//     2) static information of particle which generated the 
//        trajectory, 
//     3) trackID and parent particle ID of the trajectory,
//     4) Auxiliary points will be associated to G4SmoothTrajectoryPoint
//        to assist drawing smoothly curved trajectory.
//
// ---------------------------------------------------------------

class G4SmoothTrajectory;

#ifndef G4SmoothTrajectory_h
#define G4SmoothTrajectory_h 1

#include <stdlib.h>                   // Include from 'system'
#include <vector>

#include "trkgdefs.hh"
#include "G4VTrajectory.hh"
#include "G4Allocator.hh"
#include "G4ios.hh"                   // Include from 'system'
#include "globals.hh"                 // Include from 'global'
#include "G4ParticleDefinition.hh"    // Include from 'particle+matter'
#include "G4SmoothTrajectoryPoint.hh" // Include from 'tracking'
#include "G4Track.hh"
#include "G4Step.hh"

class G4Polyline;                     // Forward declaration

typedef std::vector<G4VTrajectoryPoint*>  TrajectoryPointContainer;

class G4SmoothTrajectory : public G4VTrajectory
{

//--------
public: // with description
//--------

// Constructor/Destrcutor

   G4SmoothTrajectory();

   G4SmoothTrajectory(const G4Track* aTrack);
   G4SmoothTrajectory(G4SmoothTrajectory &);
   virtual ~G4SmoothTrajectory();

private:
   G4SmoothTrajectory& operator= (const G4SmoothTrajectory&);

// Operators
public:
   inline void* operator new(size_t);
   inline void  operator delete(void*);
   inline int operator == (const G4SmoothTrajectory& right) const
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
   virtual void DrawTrajectory() const;
   virtual void AppendStep(const G4Step* aStep);
   virtual int GetPointEntries() const { return positionRecord->size(); }
   virtual G4VTrajectoryPoint* GetPoint(G4int i) const 
   { return (*positionRecord)[i]; }
   virtual void MergeTrajectory(G4VTrajectory* secondTrajectory);

   G4ParticleDefinition* GetParticleDefinition();

// Get method for HEPRep style attributes
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

extern G4TRACKING_DLL G4ThreadLocal
G4Allocator<G4SmoothTrajectory> *aSmoothTrajectoryAllocator;

inline void* G4SmoothTrajectory::operator new(size_t)
{
  if (!aSmoothTrajectoryAllocator)
  { aSmoothTrajectoryAllocator = new G4Allocator<G4SmoothTrajectory>; }
  return (void*)aSmoothTrajectoryAllocator->MallocSingle();
}

inline void G4SmoothTrajectory::operator delete(void* aTrajectory)
{
  aSmoothTrajectoryAllocator->FreeSingle((G4SmoothTrajectory*)aTrajectory);
}

#endif
