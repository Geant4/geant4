
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
//  The considered processes are G4MultipleScattering, G4eIonisation
//  and G4eBremsstrahlung 
//

#ifndef EMELECTRONSTANDARD_HH
#define EMELECTRONSTANDARD_HH 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"


class EMElectronStandard : public G4VPhysicsConstructor {

 public: 
   EMElectronStandard(const G4String& name = "Electron-Standard"); 
   virtual ~EMElectronStandard();
  
 protected:
   void ConstructParticle() {};
   void ConstructProcess();

   G4double facRange;
};

#endif // EMELECTRONSTANDARD_HH

