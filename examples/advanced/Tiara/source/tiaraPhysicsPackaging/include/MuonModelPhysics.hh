#ifndef MuonModelPhysics_h
#define MuonModelPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4MultipleScatteringSTD.hh"
#include "G4MuBremsstrahlungSTD.hh"
#include "G4MuPairProductionSTD.hh"
#include "G4MuIonisationSTD.hh"
#include "G4hIonisationSTD.hh"

#include "G4MuonMinusCaptureAtRest.hh"

class MuonModelPhysics : public G4VPhysicsConstructor
{
  public: 
    MuonModelPhysics(const G4String& name="muon");
    virtual ~MuonModelPhysics();

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
   G4MultipleScatteringSTD   fMuPlusMultipleScattering;
   G4MuBremsstrahlungSTD     fMuPlusBremsstrahlung ;
   G4MuPairProductionSTD     fMuPlusPairProduction;
   G4MuIonisationSTD         fMuPlusIonisation;

   G4MultipleScatteringSTD   fMuMinusMultipleScattering;
   G4MuBremsstrahlungSTD     fMuMinusBremsstrahlung ;
   G4MuPairProductionSTD     fMuMinusPairProduction;
   G4MuIonisationSTD         fMuMinusIonisation;

   G4MuonMinusCaptureAtRest fMuMinusCaptureAtRest;

   // Tau physics
   G4MultipleScatteringSTD   fTauPlusMultipleScattering;
   G4hIonisationSTD          fTauPlusIonisation;

   G4MultipleScatteringSTD   fTauMinusMultipleScattering;
   G4hIonisationSTD          fTauMinusIonisation;

};

// 2002 by J.P. Wellisch

#endif

