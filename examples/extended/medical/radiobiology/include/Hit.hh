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
/// \file radiobiology/include/Hit.hh
/// \brief Definition of the RadioBio::Hit class

#ifndef RadiobiologyHit_h
#define RadiobiologyHit_h 1

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4TouchableHandle.hh"
#include "G4VHit.hh"
#include "G4VPhysicalVolume.hh"

#include "tls.hh"

class G4ParticleDefinition;

namespace RadioBio
{

/// Detector hit class
///
/// It defines data members to store the trackID, particle type,
/// mean kinetic energy, energy deposit, initial energy, track length,
/// electronic energy, physical volume and voxel index

class Hit : public G4VHit
{
  public:
    Hit();
    Hit(const Hit&);
    ~Hit() override = default;

    // Operators
    const Hit& operator=(const Hit&);
    G4int operator==(const Hit&) const;

    inline void* operator new(size_t);
    inline void operator delete(void*);

    // Methods from base class
    void Draw() override;
    void Print() override;

    // Set methods
    void SetTrackID(G4int track) { fTrackID = track; }
    void SetPartType(const G4ParticleDefinition* part) { fParticleType = part; }
    void SetEkinMean(double EkinMean) { fEkinMean = EkinMean; }
    void SetDeltaE(double DeltaE) { fDeltaE = DeltaE; }
    void SetEinit(G4double e) { fEinit = e; }
    void SetTrackLength(G4double x) { fTrackLength = x; }
    void SetElectronEnergy(G4double elEnergy) { fElectronEnergy = elEnergy; }
    void SetPhysicalVolume(G4VPhysicalVolume* PV) { fPhysicalVolume = PV; }
    void SetVoxelIndexes(const G4TouchableHandle& TH);

    // Get methods
    G4int GetTrackID() const { return fTrackID; }
    const G4ParticleDefinition* GetPartType() const { return fParticleType; }
    G4double GetEkinMean() { return fEkinMean; }
    G4double GetDeltaE() { return fDeltaE; }
    G4double GetEinit() const { return fEinit; }
    G4double GetTrackLength() const { return fTrackLength; }
    G4double GetElectronEnergy() const { return fElectronEnergy; }
    G4VPhysicalVolume* GetPhysicalVolume() const { return fPhysicalVolume; }
    G4int GetXindex() { return fxIndex; }
    G4int GetYindex() { return fyIndex; }
    G4int GetZindex() { return fzIndex; }

  private:
    // Variables
    G4int fTrackID = -1;
    const G4ParticleDefinition* fParticleType = nullptr;
    G4double fEkinMean = 0.;
    G4double fDeltaE = 0.;
    G4double fEinit = 0.;
    G4double fTrackLength = 0.;
    G4double fElectronEnergy = 0.;
    G4VPhysicalVolume* fPhysicalVolume = nullptr;
    G4int fxIndex = -1;
    G4int fyIndex = -1;
    G4int fzIndex = -1;
};

// Definition of a HitColletion
typedef G4THitsCollection<Hit> RadioBioHitsCollection;

// Definition of the Allocator
extern G4ThreadLocal G4Allocator<Hit>* RadioBioHitAllocator;

inline void* Hit::operator new(size_t)
{
  if (!RadioBioHitAllocator) RadioBioHitAllocator = new G4Allocator<Hit>;
  return (void*)RadioBioHitAllocator->MallocSingle();
}

inline void Hit::operator delete(void* hit)
{
  RadioBioHitAllocator->FreeSingle((Hit*)hit);
}

}  // namespace RadioBio

#endif  // RadioBioHit_h
