
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
//  The considered processes are G4MultipleScattering, G4PenelopeIonisation
//  and G4PenelopeBremsstrahlung
//

#ifndef EMELECTRONPENELOPE_HH
#define EMELECTRONPENELOPE_HH 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"


class EMElectronPenelope : public G4VPhysicsConstructor {

 public: 
   EMElectronPenelope(const G4String& name = "Electron-Penelope");
   virtual ~EMElectronPenelope();
  
 protected:
   void ConstructParticle() {};
   void ConstructProcess();

   G4double facRange;
};

#endif // EMELECTRONPENELOPE_HH








