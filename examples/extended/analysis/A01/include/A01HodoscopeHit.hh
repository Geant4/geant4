// $Id: A01HodoscopeHit.hh,v 1.1 2002-11-13 07:18:31 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
#ifndef A01HodoscopeHit_h
#define A01HodoscopeHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class A01HodoscopeHit : public G4VHit
{
  public:

      A01HodoscopeHit(G4int i,G4double t);
      virtual ~A01HodoscopeHit();
      A01HodoscopeHit(const A01HodoscopeHit &right);
      const A01HodoscopeHit& operator=(const A01HodoscopeHit &right);
      int operator==(const A01HodoscopeHit &right) const;


      inline void *operator new(size_t);
      inline void operator delete(void*aHit);

      void Draw();
      void Print();

  private:
      G4int id;
      G4double time;
      G4ThreeVector pos;
      G4RotationMatrix rot;
      const G4LogicalVolume* pLogV;

  public:
      inline G4int GetID() const { return id; }
      inline G4double GetTime() const { return time; }
      inline void SetTime(G4double val) { time = val; }
      inline void SetPos(G4ThreeVector xyz) { pos = xyz; }
      inline G4ThreeVector GetPos() const { return pos; }
      inline void SetRot(G4RotationMatrix rmat) { rot = rmat; }
      inline G4RotationMatrix GetRot() const { return rot; }
      inline void SetLogV(G4LogicalVolume* val) { pLogV = val; }
      inline const G4LogicalVolume* GetLogV() const { return pLogV; }

};

typedef G4THitsCollection<A01HodoscopeHit> A01HodoscopeHitsCollection;

extern G4Allocator<A01HodoscopeHit> A01HodoscopeHitAllocator;

inline void* A01HodoscopeHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*)A01HodoscopeHitAllocator.MallocSingle();
  return aHit;
}

inline void A01HodoscopeHit::operator delete(void*aHit)
{
  A01HodoscopeHitAllocator.FreeSingle((A01HodoscopeHit*) aHit);
}

#endif


