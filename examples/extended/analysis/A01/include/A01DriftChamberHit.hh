// $Id: A01DriftChamberHit.hh,v 1.1 2002-11-13 07:17:49 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
#ifndef A01DriftChamberHit_h
#define A01DriftChamberHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class A01DriftChamberHit : public G4VHit
{
  public:

      A01DriftChamberHit();
      A01DriftChamberHit(G4int z);
      virtual ~A01DriftChamberHit();
      A01DriftChamberHit(const A01DriftChamberHit &right);
      const A01DriftChamberHit& operator=(const A01DriftChamberHit &right);
      int operator==(const A01DriftChamberHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

	  inline float x();
	  inline float y();

      virtual void Draw();
      virtual void Print();

  private:
      G4int layerID;
      G4double time;
      G4ThreeVector localPos;
      G4ThreeVector worldPos;

  public:
      inline void SetLayerID(G4int z) { layerID = z; }
      inline G4int GetLayerID() const { return layerID; }
      inline void SetTime(G4double t) { time = t; }
      inline G4double GetTime() const { return time; }
      inline void SetLocalPos(G4ThreeVector xyz) { localPos = xyz; }
      inline G4ThreeVector GetLocalPos() const { return localPos; }
      inline void SetWorldPos(G4ThreeVector xyz) { worldPos = xyz; }
      inline G4ThreeVector GetWorldPos() const { return worldPos; }
};

typedef G4THitsCollection<A01DriftChamberHit> A01DriftChamberHitsCollection;

extern G4Allocator<A01DriftChamberHit> A01DriftChamberHitAllocator;

inline void* A01DriftChamberHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*)A01DriftChamberHitAllocator.MallocSingle();
  return aHit;
}

inline void A01DriftChamberHit::operator delete(void* aHit)
{
  A01DriftChamberHitAllocator.FreeSingle((A01DriftChamberHit*) aHit);
}

#endif


