#ifndef CMSProdGeneralPhysics_h
#define CMSProdGeneralPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"


#include "G4Decay.hh"

class CMSProdGeneralPhysics : public G4VPhysicsConstructor
{
  public: 
    CMSProdGeneralPhysics(const G4String& name = "general");
    virtual ~CMSProdGeneralPhysics();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

  protected:
    G4Decay fDecayProcess;
};


#endif








