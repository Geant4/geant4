#ifndef MLLEEMPhysics_h
#define MLLEEMPhysics_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
////////////////////////////////////////////////////////////////////////////////
//
class MLLEEMPhysics : public G4VPhysicsConstructor
{
  public: 
    MLLEEMPhysics (const G4String& name ="EM");
    virtual ~MLLEEMPhysics ();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle ();

    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type
    virtual void ConstructProcess ();

};
////////////////////////////////////////////////////////////////////////////////
#endif
