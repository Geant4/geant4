#ifndef EMModelPhysics_h
#define EMModelPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4EMModelBuilder.hh"


class EMModelPhysics : public G4VPhysicsConstructor
{
  public: 
    EMModelPhysics(const G4String& name ="EM Model approach");
    virtual ~EMModelPhysics();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    G4EMModelBuilder theEMModelPhysics;
};

// 2002 by J.P. Wellisch

#endif





