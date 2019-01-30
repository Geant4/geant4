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
/// \file RE05/include/RE05CalorimeterHit.hh
/// \brief Definition of the RE05CalorimeterHit class
//

#ifndef RE05CalorimeterHit_h
#define RE05CalorimeterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;

class RE05CalorimeterHit : public G4VHit
{
  public:

      RE05CalorimeterHit();
      RE05CalorimeterHit(G4LogicalVolume* logVol,G4int z,G4int phi);
      virtual ~RE05CalorimeterHit();
      RE05CalorimeterHit(const RE05CalorimeterHit &right);
      const RE05CalorimeterHit& operator=(const RE05CalorimeterHit &right);
      G4bool operator==(const RE05CalorimeterHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      virtual void Draw();
      virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
      virtual std::vector<G4AttValue>* CreateAttValues() const;
      virtual void Print();

  private:
      G4int fZCellID;
      G4int fPhiCellID;
      G4double fEdep;
      G4ThreeVector fPos;
      G4RotationMatrix fRot;
      const G4LogicalVolume* fLogV;
      static std::map<G4String,G4AttDef> fAttDefs;

  public:
      inline void SetCellID(G4int z,G4int phi)
      {
        fZCellID = z;
        fPhiCellID = phi;
      }
      inline G4int GetZ() { return fZCellID; }
      inline G4int GetPhi() { return fPhiCellID; }
      inline void SetEdep(G4double de)
      { fEdep = de; }
      inline void AddEdep(G4double de)
      { fEdep += de; }
      inline G4double GetEdep()
      { return fEdep; }
      inline void SetPos(G4ThreeVector xyz)
      { fPos = xyz; }
      inline G4ThreeVector GetPos()
      { return fPos; }
      inline void SetRot(G4RotationMatrix rmat)
      { fRot = rmat; }
      inline G4RotationMatrix GetRot()
      { return fRot; }
      inline const G4LogicalVolume * GetLogV()
      { return fLogV; }

};

typedef G4THitsCollection<RE05CalorimeterHit> RE05CalorimeterHitsCollection;

extern G4ThreadLocal G4Allocator<RE05CalorimeterHit>* RE05CalorimeterHitAllocator;

inline void* RE05CalorimeterHit::operator new(size_t)
{
  if(!RE05CalorimeterHitAllocator)
    RE05CalorimeterHitAllocator = new G4Allocator<RE05CalorimeterHit>;
  return (void*) RE05CalorimeterHitAllocator->MallocSingle();
}

inline void RE05CalorimeterHit::operator delete(void *aHit)
{
  RE05CalorimeterHitAllocator->FreeSingle((RE05CalorimeterHit*) aHit);
}

#endif
