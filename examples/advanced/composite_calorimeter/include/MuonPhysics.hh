#ifndef MuonPhysics_h
#define MuonPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4MultipleScattering.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuIonisation.hh"
#include "G4hIonisation.hh"

#include "G4MuonMinusCaptureAtRest.hh"

class MuonPhysics : public G4VPhysicsConstructor
{
  public: 
    MuonPhysics(const G4String& name="muon");
    virtual ~MuonPhysics();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

  protected:
   // Muon physics
   G4MultipleScattering   fMuPlusMultipleScattering;
   G4MuBremsstrahlung     fMuPlusBremsstrahlung ;
   G4MuPairProduction     fMuPlusPairProduction;
   G4MuIonisation         fMuPlusIonisation;

   G4MultipleScattering   fMuMinusMultipleScattering;
   G4MuBremsstrahlung     fMuMinusBremsstrahlung ;
   G4MuPairProduction     fMuMinusPairProduction;
   G4MuIonisation         fMuMinusIonisation;

   G4MuonMinusCaptureAtRest fMuMinusCaptureAtRest;

   // Tau physics
   G4MultipleScattering   fTauPlusMultipleScattering;
   G4hIonisation          fTauPlusIonisation;

   G4MultipleScattering   fTauMinusMultipleScattering;
   G4hIonisation          fTauMinusIonisation;

};

// 2002 by J.P. Wellisch

#endif

