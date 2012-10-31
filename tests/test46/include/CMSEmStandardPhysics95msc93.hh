#ifndef SimG4Core_PhysicsLists_CMSEmStandardPhysics95msc93_h
#define SimG4Core_PhysicsLists_CMSEmStandardPhysics95msc93_h

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class CMSEmStandardPhysics95msc93 : public G4VPhysicsConstructor {

public:
  CMSEmStandardPhysics95msc93(G4int ver);
  virtual ~CMSEmStandardPhysics95msc93();

  virtual void ConstructParticle();
  virtual void ConstructProcess();

private:
  G4int               verbose;
};

#endif






