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
// $Id: A01EmCalorimeterHit.hh,v 1.3 2002-12-13 11:34:28 gunter Exp $
// --------------------------------------------------------------
//
#ifndef A01EmCalorimeterHit_h
#define A01EmCalorimeterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class A01EmCalorimeterHit : public G4VHit
{
  public:

      A01EmCalorimeterHit();
      A01EmCalorimeterHit(G4int z);
      virtual ~A01EmCalorimeterHit();
      A01EmCalorimeterHit(const A01EmCalorimeterHit &right);
      const A01EmCalorimeterHit& operator=(const A01EmCalorimeterHit &right);
      int operator==(const A01EmCalorimeterHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      virtual void Draw();
      virtual void Print();

  private:
      G4int cellID;
      G4double edep;
      G4ThreeVector pos;
      G4RotationMatrix rot;
      const G4LogicalVolume* pLogV;

  public:
      inline void SetCellID(G4int z) { cellID = z; }
      inline G4int GetCellID() const { return cellID; }
      inline void SetEdep(G4double de) { edep = de; }
      inline void AddEdep(G4double de) { edep += de; }
      inline G4double GetEdep() const { return edep; }
      inline void SetPos(G4ThreeVector xyz) { pos = xyz; }
      inline G4ThreeVector GetPos() const { return pos; }
      inline void SetRot(G4RotationMatrix rmat) { rot = rmat; }
      inline G4RotationMatrix GetRot() const { return rot; }
      inline void SetLogV(G4LogicalVolume* val) { pLogV = val; }
      inline const G4LogicalVolume* GetLogV() const { return pLogV; }
};

typedef G4THitsCollection<A01EmCalorimeterHit> A01EmCalorimeterHitsCollection;

extern G4Allocator<A01EmCalorimeterHit> A01EmCalorimeterHitAllocator;

inline void* A01EmCalorimeterHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*)A01EmCalorimeterHitAllocator.MallocSingle();
  return aHit;
}

inline void A01EmCalorimeterHit::operator delete(void* aHit)
{
  A01EmCalorimeterHitAllocator.FreeSingle((A01EmCalorimeterHit*) aHit);
}

#endif


