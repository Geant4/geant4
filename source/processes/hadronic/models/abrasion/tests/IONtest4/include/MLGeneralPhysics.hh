#ifndef MLGeneralPhysics_h
#define MLGeneralPhysics_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4Decay.hh"
////////////////////////////////////////////////////////////////////////////////
//
class MLGeneralPhysics : public G4VPhysicsConstructor
{
  public: 
    MLGeneralPhysics(const G4String& name = "general");
    virtual ~MLGeneralPhysics();

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
////////////////////////////////////////////////////////////////////////////////
//
#endif
