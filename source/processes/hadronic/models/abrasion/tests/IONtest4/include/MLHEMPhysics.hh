#ifndef MLHEMPhysics_h
#define MLHEMPhysics_h 1
/////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "G4ios.hh"
#include <vector>

#include "G4VPhysicsConstructor.hh"

#include "G4MultipleScattering.hh"

////////////////////////////////////////////////////////////////////////////////
//
class MLHEMPhysics : public G4VPhysicsConstructor
{
  public: 
    MLHEMPhysics (const G4String& name ="HEM");
    virtual ~MLHEMPhysics ();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

  private:
  
   G4String mode;

  protected:
   // Pi + 

   G4MultipleScattering thePionPlusMult;

   // Pi -
   G4MultipleScattering thePionMinusMult;

   // K + 
   G4MultipleScattering theKaonPlusMult;

   // K -
   G4MultipleScattering theKaonMinusMult;

   // Proton
   G4MultipleScattering theProtonMult;
 
   // anti-proton
   G4MultipleScattering theAntiProtonMult;
    
   // SigmaMinus 
   G4MultipleScattering theSigmaMinusMult;
  
   // AntiSigmaMinus
   G4MultipleScattering theAntiSigmaMinusMult;
   
   // SigmaPlus
   G4MultipleScattering theSigmaPlusMult;
  
   // AntiSigmaPlus
   G4MultipleScattering theAntiSigmaPlusMult;
  
   // XiMinus
   G4MultipleScattering theXiMinusMult;

   // AntiXiMinus
   G4MultipleScattering theAntiXiMinusMult;
  
   // OmegaMinus
   G4MultipleScattering theOmegaMinusMult;
   
   // AntiOmegaMinus
   G4MultipleScattering theAntiOmegaMinusMult;

};
////////////////////////////////////////////////////////////////////////////////
#endif
