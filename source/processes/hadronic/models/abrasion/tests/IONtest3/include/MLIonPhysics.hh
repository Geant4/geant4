#ifndef MLIonPhysics_h
#define MLIonPhysics_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4HadronElasticProcess.hh"
#include "G4LElastic.hh"

#include "G4DeuteronInelasticProcess.hh"
#include "G4LEDeuteronInelastic.hh"

#include "G4TritonInelasticProcess.hh"
#include "G4LETritonInelastic.hh"

#include "G4AlphaInelasticProcess.hh"
#include "G4LEAlphaInelastic.hh"

#include "G4hIonisation.hh"
#include "G4MultipleScattering.hh"
////////////////////////////////////////////////////////////////////////////////
//
class MLIonPhysics : public G4VPhysicsConstructor
{
  public: 
    MLIonPhysics (const G4String& name="Ion");
    virtual ~MLIonPhysics ();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle ();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess ();

  private:
   G4String mode;

  protected:
   // Elastic Process
   G4HadronElasticProcess       theElasticProcess;
   G4LElastic*                  theElasticModel;

   // Generic Ion physics
   G4MultipleScattering         fIonMultipleScattering;

   // Deuteron physics
   G4MultipleScattering         fDeuteronMultipleScattering;
   G4DeuteronInelasticProcess   fDeuteronProcess;
   G4LEDeuteronInelastic       *fDeuteronModel;

   // Triton physics
   G4MultipleScattering         fTritonMultipleScattering;
   G4TritonInelasticProcess     fTritonProcess;
   G4LETritonInelastic         *fTritonModel;

   // Alpha physics
   G4MultipleScattering         fAlphaMultipleScattering;
   G4AlphaInelasticProcess      fAlphaProcess;
   G4LEAlphaInelastic          *fAlphaModel;

   // He3 physics
   G4MultipleScattering         fHe3MultipleScattering;

};
////////////////////////////////////////////////////////////////////////////////
#endif
