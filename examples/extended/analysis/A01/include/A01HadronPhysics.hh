// $Id: A01HadronPhysics.hh,v 1.1 2002-11-13 07:18:24 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#ifndef A01HadronPhysics_h
#define A01HadronPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4MultipleScattering.hh"
#include "G4hIonisation.hh"

class A01HadronPhysics : public G4VPhysicsConstructor
{
  public: 
    A01HadronPhysics(const G4String& name="hadron");
    virtual ~A01HadronPhysics();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

  protected:
   // pi+/- physics
   G4MultipleScattering   fPiPlusMultipleScattering;
   G4hIonisation          fPiPlusIonisation;

   G4MultipleScattering   fPiMinusMultipleScattering;
   G4hIonisation          fPiMinusIonisation;

   // K+/- physics
   G4MultipleScattering   fKaonPlusMultipleScattering;
   G4hIonisation          fKaonPlusIonisation;

   G4MultipleScattering   fKaonMinusMultipleScattering;
   G4hIonisation          fKaonMinusIonisation;

};


#endif

