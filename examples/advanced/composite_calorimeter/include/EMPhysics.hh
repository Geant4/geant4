#ifndef EMPhysics_h
#define EMPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4EMBuilder.hh"


class EMPhysics : public G4VPhysicsConstructor
{
  public: 
    EMPhysics(const G4String& name ="EM");
    virtual ~EMPhysics();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    G4EMBuilder theEMPhysics;
};

// 2002 by J.P. Wellisch

#endif





