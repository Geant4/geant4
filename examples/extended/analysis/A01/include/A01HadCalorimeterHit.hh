// $Id: A01HadCalorimeterHit.hh,v 1.1 2002-11-13 07:18:21 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
#ifndef A01HadCalorimeterHit_h
#define A01HadCalorimeterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class A01HadCalorimeterHit : public G4VHit
{
  public:

      A01HadCalorimeterHit();
      A01HadCalorimeterHit(G4int iCol,G4int iRow);
      virtual ~A01HadCalorimeterHit();
      A01HadCalorimeterHit(const A01HadCalorimeterHit &right);
      const A01HadCalorimeterHit& operator=(const A01HadCalorimeterHit &right);
      int operator==(const A01HadCalorimeterHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      virtual void Draw();
      virtual void Print();

  private:
      G4int columnID;
      G4int rowID;
      G4double edep;
      G4ThreeVector pos;
      G4RotationMatrix rot;

  public:
      inline void SetColumnID(G4int z) { columnID = z; }
      inline G4int GetColumnID() const { return columnID; }
      inline void SetRowID(G4int z) { rowID = z; }
      inline G4int GetRowID() const { return rowID; }
      inline void SetEdep(G4double de) { edep = de; }
      inline void AddEdep(G4double de) { edep += de; }
      inline G4double GetEdep() const { return edep; }
      inline void SetPos(G4ThreeVector xyz) { pos = xyz; }
      inline G4ThreeVector GetPos() const { return pos; }
      inline void SetRot(G4RotationMatrix rmat) { rot = rmat; }
      inline G4RotationMatrix GetRot() const { return rot; }

};

typedef G4THitsCollection<A01HadCalorimeterHit> A01HadCalorimeterHitsCollection;

extern G4Allocator<A01HadCalorimeterHit> A01HadCalorimeterHitAllocator;

inline void* A01HadCalorimeterHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*)A01HadCalorimeterHitAllocator.MallocSingle();
  return aHit;
}

inline void A01HadCalorimeterHit::operator delete(void* aHit)
{
  A01HadCalorimeterHitAllocator.FreeSingle((A01HadCalorimeterHit*) aHit);
}

#endif


