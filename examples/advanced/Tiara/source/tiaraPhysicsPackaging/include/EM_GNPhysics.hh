#ifndef EM_GNPhysics_h
#define EM_GNPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4EMBuilder.hh"
#include "G4ElectroNuclearBuilder.hh"


class EM_GNPhysics : public G4VPhysicsConstructor
{
  public: 
    EM_GNPhysics(const G4String& name ="EM");
    virtual ~EM_GNPhysics();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    G4EMBuilder theEMPhysics;
    G4ElectroNuclearBuilder theGNPhysics;
};

// 2002 by J.P. Wellisch

#endif





