//========================
// taken from example N04 of G4 distribution
//========================
#ifndef IonModelPhysics_h
#define IonModelPhysics_h 1

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

#include "G4hIonisationSTD.hh"
#include "G4ionIonisation.hh"
#include "G4MultipleScatteringSTD.hh"

class IonModelPhysics : public G4VPhysicsConstructor
{
  public: 
    IonModelPhysics(const G4String& name="ion");
    virtual ~IonModelPhysics();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

  protected:
   // Elastic Process
   G4LElastic*            theElasticModel;

   // Generic Ion physics
   G4HadronElasticProcess    theIonElasticProcess;
   G4MultipleScatteringSTD   fIonMultipleScattering;
   G4ionIonisation           fIonIonisation;

   // Deuteron physics
   G4HadronElasticProcess      theDElasticProcess;
   G4MultipleScatteringSTD     fDeuteronMultipleScattering;
   G4hIonisationSTD            fDeuteronIonisation;
   G4DeuteronInelasticProcess  fDeuteronProcess;
   G4LEDeuteronInelastic*      fDeuteronModel;

   // Triton physics
   G4HadronElasticProcess      theTElasticProcess;
   G4MultipleScatteringSTD     fTritonMultipleScattering;
   G4hIonisationSTD            fTritonIonisation;
   G4TritonInelasticProcess    fTritonProcess;
   G4LETritonInelastic*        fTritonModel;
  
   // Alpha physics
   G4HadronElasticProcess      theAElasticProcess;
   G4MultipleScatteringSTD     fAlphaMultipleScattering;
   G4hIonisationSTD            fAlphaIonisation;
   G4AlphaInelasticProcess     fAlphaProcess;
   G4LEAlphaInelastic*         fAlphaModel;

   // He3 physics
   G4HadronElasticProcess      theHe3ElasticProcess;
   G4MultipleScatteringSTD     fHe3MultipleScattering;
   G4hIonisationSTD            fHe3Ionisation;

};

// 2002 by J.P. Wellisch

#endif

