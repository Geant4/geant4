// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN04IonPhysics.hh,v 1.1 2001-01-06 06:56:16 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
// Class Description:
//      This class is an derived class of G4VPhysicsConstructor
//
// ------------------------------------------- 
//	History
//        first version                   12 Nov. 2000 by H.Kurashige 
// ------------------------------------------------------------
#ifndef ExN04IonPhysics_h
#define ExN04IonPhysics_h 1

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

class ExN04IonPhysics : public G4VPhysicsConstructor
{
  public: 
    ExN04IonPhysics(const G4String& name="ion");
    virtual ~ExN04IonPhysics();

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
   G4HadronElasticProcess theElasticProcess;
   G4LElastic*            theElasticModel;

   // Generic Ion physics
   G4MultipleScattering   fIonMultipleScattering;
   G4hIonisation          fIonIonisation;

   // Deuteron physics
   G4MultipleScattering        fDeuteronMultipleScattering;
   G4hIonisation               fDeuteronIonisation;
   G4DeuteronInelasticProcess  fDeuteronProcess;
   G4LEDeuteronInelastic*      fDeuteronModel;

   // Triton physics
   G4MultipleScattering        fTritonMultipleScattering;
   G4hIonisation               fTritonIonisation;
   G4TritonInelasticProcess    fTritonProcess;
   G4LETritonInelastic*        fTritonModel;
  
   // Alpha physics
   G4MultipleScattering        fAlphaMultipleScattering;
   G4hIonisation               fAlphaIonisation;
   G4AlphaInelasticProcess     fAlphaProcess;
   G4LEAlphaInelastic*         fAlphaModel;

   // He3 physics
   G4MultipleScattering        fHe3MultipleScattering;
   G4hIonisation               fHe3Ionisation;

};


#endif

