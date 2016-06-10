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
// $Id: RE05CalorimeterHit.hh 69764 2013-05-14 09:59:36Z gcosmo $
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
      G4int operator==(const RE05CalorimeterHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      virtual void Draw();
      virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
      virtual std::vector<G4AttValue>* CreateAttValues() const;
      virtual void Print();

  private:
      G4int ZCellID;
      G4int PhiCellID;
      G4double edep;
      G4ThreeVector pos;
      G4RotationMatrix rot;
      const G4LogicalVolume* pLogV;
      static std::map<G4String,G4AttDef> fAttDefs;

  public:
      inline void SetCellID(G4int z,G4int phi)
      {
        ZCellID = z;
        PhiCellID = phi;
      }
      inline G4int GetZ() { return ZCellID; }
      inline G4int GetPhi() { return PhiCellID; }
      inline void SetEdep(G4double de)
      { edep = de; }
      inline void AddEdep(G4double de)
      { edep += de; }
      inline G4double GetEdep()
      { return edep; }
      inline void SetPos(G4ThreeVector xyz)
      { pos = xyz; }
      inline G4ThreeVector GetPos()
      { return pos; }
      inline void SetRot(G4RotationMatrix rmat)
      { rot = rmat; }
      inline G4RotationMatrix GetRot()
      { return rot; }
      inline const G4LogicalVolume * GetLogV()
      { return pLogV; }

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


