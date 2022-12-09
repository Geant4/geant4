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
/// \file B5/include/EmCalorimeterHit.hh
/// \brief Definition of the B5::EmCalorimeterHit class

#ifndef B5EmCalorimeterHit_h
#define B5EmCalorimeterHit_h 1

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

/// EM Calorimeter hit
///
/// It records:
/// - the cell ID
/// - the energy deposit
/// - the cell logical volume, its position and rotation

class EmCalorimeterHit : public G4VHit
{
  public:
    EmCalorimeterHit() = default;
    EmCalorimeterHit(G4int cellID);
    EmCalorimeterHit(const EmCalorimeterHit &right) = default;
    ~EmCalorimeterHit() override = default;

    EmCalorimeterHit& operator=(const EmCalorimeterHit &right) = default;
    G4bool operator==(const EmCalorimeterHit &right) const;

    inline void *operator new(size_t);
    inline void operator delete(void *aHit);

    void Draw() override;
    const std::map<G4String,G4AttDef>* GetAttDefs() const override;
    std::vector<G4AttValue>* CreateAttValues() const override;
    void Print() override;

    void SetCellID(G4int z) { fCellID = z; }
    G4int GetCellID() const { return fCellID; }

    void SetEdep(G4double de) { fEdep = de; }
    void AddEdep(G4double de) { fEdep += de; }
    G4double GetEdep() const { return fEdep; }

    void SetPos(G4ThreeVector xyz) { fPos = xyz; }
    G4ThreeVector GetPos() const { return fPos; }

    void SetRot(G4RotationMatrix rmat) { fRot = rmat; }
    G4RotationMatrix GetRot() const { return fRot; }

    void SetLogV(G4LogicalVolume* val) { fPLogV = val; }
    const G4LogicalVolume* GetLogV() const { return fPLogV; }

  private:
    G4int fCellID = -1;
    G4double fEdep = 0.;
    G4ThreeVector fPos;
    G4RotationMatrix fRot;
    const G4LogicalVolume* fPLogV = nullptr;
};

using EmCalorimeterHitsCollection = G4THitsCollection<EmCalorimeterHit>;

extern G4ThreadLocal G4Allocator<EmCalorimeterHit>* EmCalorimeterHitAllocator;

inline void* EmCalorimeterHit::operator new(size_t)
{
  if (!EmCalorimeterHitAllocator) {
       EmCalorimeterHitAllocator = new G4Allocator<EmCalorimeterHit>;
  }
  return (void*)EmCalorimeterHitAllocator->MallocSingle();
}

inline void EmCalorimeterHit::operator delete(void* aHit)
{
  EmCalorimeterHitAllocator->FreeSingle((EmCalorimeterHit*) aHit);
}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
