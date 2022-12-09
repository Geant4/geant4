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
/// \file B5/include/HadCalorimeterHit.hh
/// \brief Definition of the B5::HadCalorimeterHit class

#ifndef B5HadCalorimeterHit_h
#define B5HadCalorimeterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

namespace B5
{

/// Hadron Calorimeter hit
///
/// It records:
/// - the cell column ID and row ID
/// - the energy deposit
/// - the cell position and rotation

class HadCalorimeterHit : public G4VHit
{
  public:
    HadCalorimeterHit() = default;
    HadCalorimeterHit(G4int iCol,G4int iRow);
    HadCalorimeterHit(const HadCalorimeterHit &right) = default;
    ~HadCalorimeterHit() override = default;

    HadCalorimeterHit& operator=(const HadCalorimeterHit &right) = default;
    G4bool operator==(const HadCalorimeterHit &right) const;

    inline void *operator new(size_t);
    inline void operator delete(void *aHit);

    void Draw() override;
    const std::map<G4String,G4AttDef>* GetAttDefs() const override;
    std::vector<G4AttValue>* CreateAttValues() const override;
    void Print() override;

    void SetColumnID(G4int z) { fColumnID = z; }
    G4int GetColumnID() const { return fColumnID; }

    void SetRowID(G4int z) { fRowID = z; }
    G4int GetRowID() const { return fRowID; }

    void SetEdep(G4double de) { fEdep = de; }
    void AddEdep(G4double de) { fEdep += de; }
    G4double GetEdep() const { return fEdep; }

    void SetPos(G4ThreeVector xyz) { fPos = xyz; }
    G4ThreeVector GetPos() const { return fPos; }

    void SetRot(G4RotationMatrix rmat) { fRot = rmat; }
    G4RotationMatrix GetRot() const { return fRot; }

  private:
    G4int fColumnID = -1;
    G4int fRowID = -1;
    G4double fEdep = 0.;
    G4ThreeVector fPos;
    G4RotationMatrix fRot;
};

using HadCalorimeterHitsCollection = G4THitsCollection<HadCalorimeterHit>;

extern G4ThreadLocal G4Allocator<HadCalorimeterHit>* HadCalorimeterHitAllocator;

inline void* HadCalorimeterHit::operator new(size_t)
{
  if (!HadCalorimeterHitAllocator) {
       HadCalorimeterHitAllocator = new G4Allocator<HadCalorimeterHit>;
  }
  return (void*)HadCalorimeterHitAllocator->MallocSingle();
}

inline void HadCalorimeterHit::operator delete(void* aHit)
{
  HadCalorimeterHitAllocator->FreeSingle((HadCalorimeterHit*) aHit);
}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
