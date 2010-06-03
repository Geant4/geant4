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
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef F04Trajectory_h
#define F04Trajectory_h 1

#include <vector>
#include <stdlib.h>

#include "globals.hh"

#include "G4ios.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VTrajectory.hh"
#include "G4ParticleDefinition.hh" 
#include "G4TrajectoryPoint.hh"

typedef std::vector<G4VTrajectoryPoint*> TrajectoryPointContainer;

class F04Trajectory : public G4VTrajectory
{

//--------
   public: // without description
//--------

// Constructor/Destructor

     F04Trajectory();
     F04Trajectory(const G4Track* aTrack);
     F04Trajectory(F04Trajectory &);
     virtual ~F04Trajectory();

// Operators

     inline void* operator new(size_t);
     inline void  operator delete(void*);
     inline int operator == (const F04Trajectory& right) const
     { return (this==&right); }

// Get/Set functions 

     inline G4int GetTrackID() const { return fTrackID; }
     inline G4int GetParentID() const { return fParentID; }
     inline G4String GetParticleName() const { return ParticleName; }
     inline G4double GetCharge() const { return PDGCharge; }
     inline G4int GetPDGEncoding() const { return PDGEncoding; }
     inline G4ThreeVector GetInitialMomentum() const { return initialMomentum; }

// Other member functions

     virtual void ShowTrajectory(std::ostream& os=G4cout) const;
     virtual void DrawTrajectory() const;
     virtual void DrawTrajectory(G4int i_mode=0) const;
     virtual void AppendStep(const G4Step* aStep);
     virtual void MergeTrajectory(G4VTrajectory* secondTrajectory);

     G4ParticleDefinition* GetParticleDefinition();

     virtual int GetPointEntries() const 
     { return fpPointsContainer->size(); }
     virtual G4VTrajectoryPoint* GetPoint(G4int i) const
     { return (*fpPointsContainer)[i]; }

    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;

//---------
   private:
//---------

// Member data

     TrajectoryPointContainer* fpPointsContainer;

     G4int fTrackID;
     G4int fParentID;
     G4double PDGCharge;
     G4int    PDGEncoding;
     G4String ParticleName;
     G4ThreeVector initialMomentum;

};

#if defined G4TRACKING_ALLOC_EXPORT
  extern G4DLLEXPORT G4Allocator<F04Trajectory> aTrajectoryAllocator;
#else
  extern G4DLLIMPORT G4Allocator<F04Trajectory> aTrajectoryAllocator;
#endif

inline void* F04Trajectory::operator new(size_t) {
    void* aTrajectory = (void*) aTrajectoryAllocator.MallocSingle();
    return aTrajectory;
}

inline void F04Trajectory::operator delete(void* aTrajectory) {
    aTrajectoryAllocator.FreeSingle((F04Trajectory*)aTrajectory);
}

#endif
