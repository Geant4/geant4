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
/// \file B5DriftChamberHit.hh
/// \brief Definition of the B5DriftChamberHit class

#ifndef B5DriftChamberHit_h
#define B5DriftChamberHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

/// Drift chamber hit
///
/// It records:
/// - the layer ID
/// - the particle time
/// - the particle local and global positions

class B5DriftChamberHit : public G4VHit
{
  public:
    B5DriftChamberHit();
    B5DriftChamberHit(G4int layerID);
    B5DriftChamberHit(const B5DriftChamberHit &right);
    virtual ~B5DriftChamberHit();

    const B5DriftChamberHit& operator=(const B5DriftChamberHit &right);
    G4bool operator==(const B5DriftChamberHit &right) const;
    
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
    
    virtual void Draw();
    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;
    virtual void Print();

    void SetLayerID(G4int z) { fLayerID = z; }
    G4int GetLayerID() const { return fLayerID; }

    void SetTime(G4double t) { fTime = t; }
    G4double GetTime() const { return fTime; }

    void SetLocalPos(G4ThreeVector xyz) { fLocalPos = xyz; }
    G4ThreeVector GetLocalPos() const { return fLocalPos; }

    void SetWorldPos(G4ThreeVector xyz) { fWorldPos = xyz; }
    G4ThreeVector GetWorldPos() const { return fWorldPos; }
    
  private:
    G4int fLayerID;
    G4double fTime;
    G4ThreeVector fLocalPos;
    G4ThreeVector fWorldPos;
};

using B5DriftChamberHitsCollection = G4THitsCollection<B5DriftChamberHit>;

extern G4ThreadLocal G4Allocator<B5DriftChamberHit>* B5DriftChamberHitAllocator;

inline void* B5DriftChamberHit::operator new(size_t)
{
  if (!B5DriftChamberHitAllocator) {
       B5DriftChamberHitAllocator = new G4Allocator<B5DriftChamberHit>;
  }
  return (void*)B5DriftChamberHitAllocator->MallocSingle();
}

inline void B5DriftChamberHit::operator delete(void* aHit)
{
  B5DriftChamberHitAllocator->FreeSingle((B5DriftChamberHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
