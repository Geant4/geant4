#ifndef SimG4Core_PhysicsLists_CMSEmStandardPhysics92_h
#define SimG4Core_PhysicsLists_CMSEmStandardPhysics92_h

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class CMSEmStandardPhysics92 : public G4VPhysicsConstructor {

public:
  CMSEmStandardPhysics92(G4int ver);
  virtual ~CMSEmStandardPhysics92();

  virtual void ConstructParticle();
  virtual void ConstructProcess();

private:
  G4int               verbose;
};

#endif






