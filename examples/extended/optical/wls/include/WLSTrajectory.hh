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
/// \file optical/wls/include/WLSTrajectory.hh
/// \brief Definition of the WLSTrajectory class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef WLSTrajectory_h
#define WLSTrajectory_h 1

#include "globals.hh"
#include "G4Allocator.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4TrajectoryPoint.hh"
#include "G4VTrajectory.hh"

#include <vector>

typedef std::vector<G4VTrajectoryPoint*> WLSTrajectoryPointContainer;

class WLSTrajectory : public G4VTrajectory
{
 public:
  WLSTrajectory() = default;
  WLSTrajectory(const G4Track*);
  WLSTrajectory(WLSTrajectory&);
  ~WLSTrajectory() override;

  inline void* operator new(size_t);
  inline void operator delete(void*);
  inline int operator==(const WLSTrajectory& right) const
  {
    return (this == &right);
  }

  inline G4int GetTrackID() const override { return fTrackID; }
  inline G4int GetParentID() const override { return fParentID; }
  inline G4String GetParticleName() const override { return fParticleName; }
  inline G4double GetCharge() const override { return fPDGCharge; }
  inline G4int GetPDGEncoding() const override { return fPDGEncoding; }
  inline G4ThreeVector GetInitialMomentum() const override
  {
    return fInitialMomentum;
  }

  void ShowTrajectory(std::ostream& os = G4cout) const override;
  void AppendStep(const G4Step* aStep) override;
  void MergeTrajectory(G4VTrajectory* secondTrajectory) override;

  G4ParticleDefinition* GetParticleDefinition();

  int GetPointEntries() const override { return fpPointsContainer->size(); }
  G4VTrajectoryPoint* GetPoint(G4int i) const override
  {
    return (*fpPointsContainer)[i];
  }

  const std::map<G4String, G4AttDef>* GetAttDefs() const override;
  std::vector<G4AttValue>* CreateAttValues() const override;

 private:
  WLSTrajectoryPointContainer* fpPointsContainer = nullptr;

  G4int fTrackID = 0;
  G4int fParentID = 0;
  G4double fPDGCharge = 0.;
  G4int fPDGEncoding = 0;
  G4String fParticleName;
  G4ThreeVector fInitialMomentum;

  G4ParticleDefinition* fParticleDefinition = nullptr;
};

extern G4ThreadLocal G4Allocator<WLSTrajectory>* WLSTrajectoryAllocator;

inline void* WLSTrajectory::operator new(size_t)
{
  if(!WLSTrajectoryAllocator)
    WLSTrajectoryAllocator = new G4Allocator<WLSTrajectory>;
  return (void*) WLSTrajectoryAllocator->MallocSingle();
}

inline void WLSTrajectory::operator delete(void* aTrajectory)
{
  WLSTrajectoryAllocator->FreeSingle((WLSTrajectory*) aTrajectory);
}

#endif
