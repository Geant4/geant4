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
/// \file field/field04/include/F04Trajectory.hh
/// \brief Definition of the F04Trajectory class
//

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

using TrajectoryPointContainer = std::vector<G4VTrajectoryPoint *>;

class F04Trajectory : public G4VTrajectory
{

//--------
   public: // without description
//--------

// Constructor/Destructor

     F04Trajectory();
     F04Trajectory(const G4Track* aTrack);
     F04Trajectory(F04Trajectory&);
     ~F04Trajectory() override;

// Operators

     inline void* operator new(size_t);
     inline void  operator delete(void*);
     inline int operator == (const F04Trajectory& right) const
     { return (this==&right); }

// Get/Set functions

     inline G4int GetTrackID() const override { return fTrackID; }
     inline G4int GetParentID() const override { return fParentID; }
     inline G4String GetParticleName() const override { return fParticleName; }
     inline G4double GetCharge() const override { return fPDGCharge; }
     inline G4int GetPDGEncoding() const override { return fPDGEncoding; }
     inline G4ThreeVector GetInitialMomentum() const override {return fInitialMomentum;}

// Other member functions

     void ShowTrajectory(std::ostream& os=G4cout) const override;
     void AppendStep(const G4Step* aStep) override;
     void MergeTrajectory(G4VTrajectory* secondTrajectory) override;

     G4ParticleDefinition* GetParticleDefinition();

     int GetPointEntries() const override
     { return fpPointsContainer->size(); }
     G4VTrajectoryPoint* GetPoint(G4int i) const override
     { return (*fpPointsContainer)[i]; }

    const std::map<G4String,G4AttDef>* GetAttDefs() const override;
    std::vector<G4AttValue>* CreateAttValues() const override;

//---------
   private:
//---------

// Member data

     TrajectoryPointContainer* fpPointsContainer = nullptr;

     G4int         fTrackID = 0;
     G4int         fParentID = 0;
     G4double      fPDGCharge = 0.;
     G4int         fPDGEncoding = 0;
     G4String      fParticleName;
     G4ThreeVector fInitialMomentum;

};

extern G4ThreadLocal G4Allocator<F04Trajectory>* F04TrajectoryAllocator;

inline void* F04Trajectory::operator new(size_t) {
    if(!F04TrajectoryAllocator)
      F04TrajectoryAllocator = new G4Allocator<F04Trajectory>;
    return (void*) F04TrajectoryAllocator->MallocSingle();
}

inline void F04Trajectory::operator delete(void* aTrajectory) {
    F04TrajectoryAllocator->FreeSingle((F04Trajectory*)aTrajectory);
}

#endif
