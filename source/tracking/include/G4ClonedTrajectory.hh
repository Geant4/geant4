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
// G4ClonedTrajectory
//
// Class description:
//
// This class is identical to G4Trajectory class, but uses the
// singleton G4Allocator so that the master thread can safely
// delete the object instantiated by worker threads in sub-event
// parallel mode. This class object should be instantiated as
// a clone of G4Trajectory object.
//
// Makoto Asai (JLab) - Oct.2024
// --------------------------------------------------------------------
//
#ifndef G4ClonedTrajectory_hh
#define G4ClonedTrajectory_hh 1

#include "G4Allocator.hh"
#include "G4ParticleDefinition.hh"  // Include from 'particle+matter'
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ClonedTrajectoryPoint.hh"  // Include from 'tracking'
#include "G4VTrajectory.hh"
#include "G4ios.hh"  // Include from 'system'
#include "globals.hh"  // Include from 'global'
#include "G4Threading.hh"

#include "trkgdefs.hh"
#include <stdlib.h>  // Include from 'system'

#include <vector>

class G4Polyline;
class G4Trajectory;

class G4ClonedTrajectory : public G4VTrajectory
{
  using G4ClonedTrajectoryPointContainer = std::vector<G4VTrajectoryPoint*>;

 public:
  // Constructors/Destructor

  G4ClonedTrajectory() = default;
  G4ClonedTrajectory(const G4Trajectory&);
  ~G4ClonedTrajectory() override;

  // Operators

  inline void* operator new(size_t);
  inline void operator delete(void*);
  inline G4bool operator==(const G4ClonedTrajectory& r) const;

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
  G4ClonedTrajectoryPointContainer* positionRecord = nullptr;
  G4int fTrackID = 0;
  G4int fParentID = 0;
  G4int PDGEncoding = 0;
  G4double PDGCharge = 0.0;
  G4String ParticleName = "dummy";
  G4double initialKineticEnergy = 0.0;
  G4ThreeVector initialMomentum;
};

extern G4TRACKING_DLL G4Allocator<G4ClonedTrajectory>*& aClonedTrajectoryAllocator();

inline void* G4ClonedTrajectory::operator new(size_t)
{
  if (aClonedTrajectoryAllocator() == nullptr) {
    aClonedTrajectoryAllocator() = new G4Allocator<G4ClonedTrajectory>;
  }
  return (void*)aClonedTrajectoryAllocator()->MallocSingle();
}

inline void G4ClonedTrajectory::operator delete(void* aTrajectory)
{
  aClonedTrajectoryAllocator()->FreeSingle((G4ClonedTrajectory*)aTrajectory);
}

inline G4bool G4ClonedTrajectory::operator==(const G4ClonedTrajectory& r) const { return (this == &r); }

#endif
