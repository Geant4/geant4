
// ==============================
// PHYSICS PROCESSES:
// ==============================
//  Electromagnetic processes for: 
//      electrons
//     
// 
// ==============================
// COMMENTS:
// ==============================
//  The considered processes are G4eMultipleScattering, G4LivermoreIonisation
//  and G4LivemoreBremsstrahlung (Low Energy Package, EEDL libraries)
//

#ifndef EMELECTRONEEDL_HH
#define EMELECTRONEEDL_HH 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"


class EMElectronEEDL : public G4VPhysicsConstructor {

 public: 
   EMElectronEEDL(const G4String& name = "Electron-EEDL"); 
   virtual ~EMElectronEEDL();
  
 protected:
   void ConstructParticle() {}; 
   void ConstructProcess();

   G4double facRange;
};

#endif // EMELECTRONEEDL_HH








