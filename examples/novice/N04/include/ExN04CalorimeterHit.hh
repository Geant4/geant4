//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//

#ifndef ExN04CalorimeterHit_h
#define ExN04CalorimeterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class ExN04CalorimeterHit : public G4VHit
{
  public:

      ExN04CalorimeterHit();
      ExN04CalorimeterHit(G4LogicalVolume* logVol,G4int z,G4int phi);
      ~ExN04CalorimeterHit();
      ExN04CalorimeterHit(const ExN04CalorimeterHit &right);
      const ExN04CalorimeterHit& operator=(const ExN04CalorimeterHit &right);
      int operator==(const ExN04CalorimeterHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      void Draw();
      void Print();

  private:
      G4int ZCellID;
      G4int PhiCellID;
      G4double edep;
      G4ThreeVector pos;
      G4RotationMatrix rot;
      const G4LogicalVolume* pLogV;

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

typedef G4THitsCollection<ExN04CalorimeterHit> ExN04CalorimeterHitsCollection;

extern G4Allocator<ExN04CalorimeterHit> ExN04CalorimeterHitAllocator;

inline void* ExN04CalorimeterHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) ExN04CalorimeterHitAllocator.MallocSingle();
  return aHit;
}

inline void ExN04CalorimeterHit::operator delete(void *aHit)
{
  ExN04CalorimeterHitAllocator.FreeSingle((ExN04CalorimeterHit*) aHit);
}

#endif


