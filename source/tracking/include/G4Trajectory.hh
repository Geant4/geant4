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
// G4Trajectory
//
// Class description:
//
// This class represents the trajectory of a particle being tracked.
// It includes information of:
//     1) List of trajectory points which compose the trajectory;
//     2) Static information of the particle which generated the
//        trajectory;
//     3) Track ID and parent particle ID of the trajectory.

// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Makoto  Asai   (e-mail: asai@slac.stanford.edu)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
// --------------------------------------------------------------------
#ifndef G4Trajectory_hh
#define G4Trajectory_hh 1

#include "G4Allocator.hh"
#include "G4ParticleDefinition.hh"  // Include from 'particle+matter'
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4TrajectoryPoint.hh"  // Include from 'tracking'
#include "G4VTrajectory.hh"
#include "G4ios.hh"  // Include from 'system'
#include "globals.hh"  // Include from 'global'

#include "trkgdefs.hh"
#include <stdlib.h>  // Include from 'system'

#include <vector>

class G4Polyline;

class G4Trajectory : public G4VTrajectory
{
  using G4TrajectoryPointContainer = std::vector<G4VTrajectoryPoint*>;

 public:
  // Constructors/Destructor

  G4Trajectory() = default;
  G4Trajectory(const G4Track* aTrack);
  G4Trajectory(G4Trajectory&);
  ~G4Trajectory() override;

  // Operators

  inline void* operator new(size_t);
  inline void operator delete(void*);
  inline G4bool operator==(const G4Trajectory& r) const;

  // Get/Set functions

  inline G4int GetTrackID() const override { return fTrackID; }
  inline G4int GetParentID() const override { return fParentID; }
  inline G4String GetParticleName() const override { return ParticleName; }
  inline G4double GetCharge() const override { return PDGCharge; }
  inline G4int GetPDGEncoding() const override { return PDGEncoding; }
  inline G4double GetInitialKineticEnergy() const { return initialKineticEnergy; }
  inline G4ThreeVector GetInitialMomentum() const override { return initialMomentum; }

  // Other member functions

  void ShowTrajectory(std::ostream& os = G4cout) const override;
  void DrawTrajectory() const override;
  void AppendStep(const G4Step* aStep) override;
  G4int GetPointEntries() const override { return G4int(positionRecord->size()); }
  G4VTrajectoryPoint* GetPoint(G4int i) const override { return (*positionRecord)[i]; }
  void MergeTrajectory(G4VTrajectory* secondTrajectory) override;

  G4ParticleDefinition* GetParticleDefinition();

  const std::map<G4String, G4AttDef>* GetAttDefs() const override;
  std::vector<G4AttValue>* CreateAttValues() const override;

 private:
  G4TrajectoryPointContainer* positionRecord = nullptr;
  G4int fTrackID = 0;
  G4int fParentID = 0;
  G4int PDGEncoding = 0;
  G4double PDGCharge = 0.0;
  G4String ParticleName = "";
  G4double initialKineticEnergy = 0.0;
  G4ThreeVector initialMomentum;
};

extern G4TRACKING_DLL G4Allocator<G4Trajectory>*& aTrajectoryAllocator();

inline void* G4Trajectory::operator new(size_t)
{
  if (aTrajectoryAllocator() == nullptr) {
    aTrajectoryAllocator() = new G4Allocator<G4Trajectory>;
  }
  return (void*)aTrajectoryAllocator()->MallocSingle();
}

inline void G4Trajectory::operator delete(void* aTrajectory)
{
  aTrajectoryAllocator()->FreeSingle((G4Trajectory*)aTrajectory);
}

inline G4bool G4Trajectory::operator==(const G4Trajectory& r) const { return (this == &r); }

#endif
