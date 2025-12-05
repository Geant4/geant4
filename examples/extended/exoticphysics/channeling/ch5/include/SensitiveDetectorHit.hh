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
// gpaterno, October 2025
//
/// \file SensitiveDetectorHit.hh
/// \brief Definition of the SensitiveDetectorHit class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef SensitiveDetectorHit_h
#define SensitiveDetectorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// SensitiveDetectorHit class.

class G4AttDef;
class G4AttValue;

class SensitiveDetectorHit : public G4VHit
{
public:
    SensitiveDetectorHit() = default;
    ~SensitiveDetectorHit() override = default;
    
    SensitiveDetectorHit(const SensitiveDetectorHit &right) = default;
    SensitiveDetectorHit& operator=(const SensitiveDetectorHit &right) = default;
    
    int operator==(const SensitiveDetectorHit &right) const;
    
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
    
    void Draw() override;
    const std::map<G4String,G4AttDef>* GetAttDefs() const override;
    std::vector<G4AttValue>* CreateAttValues() const override;
    void Print() override;

    inline void SetTrackID(G4int z) {fTrackID = z;}
    inline G4int GetTrackID() const {return fTrackID;}

    inline void SetTrackIDP(G4int z) {fTrackIDP = z;}
    inline G4int GetTrackIDP() const {return fTrackIDP;}

    inline void SetTime(G4double t) {fTime = t;}
    inline G4double GetTime() const {return fTime;}
    
    inline void SetPos(G4ThreeVector xyz) {fPos = xyz;}
    inline G4ThreeVector GetPos() const {return fPos;}
    
    inline void SetMom(G4ThreeVector xyz) {fMom = xyz;}
    inline G4ThreeVector GetMom() const {return fMom;}
    
    inline void SetEnergy(G4double energy) {fEnergy = energy;}
    inline G4double GetEnergy() const {return fEnergy;}
       
    inline void SetType(G4int type) {fType = type;}
    inline G4int GetType() const {return fType;}  
      
    inline void SetParticle(G4String particle) {fParticle = particle;}
    inline G4String GetParticle() const {return fParticle;}

    inline void SetWeight(G4double weight) {fWeight = weight;}
    inline G4double GetWeight() const {return fWeight;}

    inline void SetDetID(G4int did) {fDetID = did;}
    inline G4int GetDetID() const {return fDetID;}
    
private:
    G4int fTrackID = -1;
    G4int fTrackIDP = -1;
    G4double fTime = 0;
    G4ThreeVector fPos = G4ThreeVector(0., 0., 0.);
    G4ThreeVector fMom = G4ThreeVector(0., 0., 0.);
    G4double fEnergy = 0.;
    G4int fType = -11;
    G4String fParticle = "";
    G4double fWeight = -1;
    G4int fDetID = -1;
};

using SensitiveDetectorHitsCollection = G4THitsCollection<SensitiveDetectorHit>;

extern G4ThreadLocal G4Allocator<SensitiveDetectorHit>* SensitiveDetectorHitAllocator;

inline void* SensitiveDetectorHit::operator new(size_t)
{
  if (!SensitiveDetectorHitAllocator) SensitiveDetectorHitAllocator = 
      new G4Allocator<SensitiveDetectorHit>;
  return (void *) SensitiveDetectorHitAllocator->MallocSingle();
}

inline void SensitiveDetectorHit::operator delete(void* aHit)
{
  SensitiveDetectorHitAllocator->FreeSingle((SensitiveDetectorHit*) aHit);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
