
// ==============================
// PHYSICS PROCESSES:
// ==============================
//  Electromagnetic processes for: 
//      photons
//     
// 
// ==============================
// COMMENTS:
// ==============================
//  The considered processes are G4PhotoElectricEffect, G4ComptonScattering
//  and G4GammaConversion 
//

#ifndef EMPHOTONSTANDARD_HH
#define EMPHOTONSTANDARD_HH 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"


class EMPhotonStandard : public G4VPhysicsConstructor {

 public: 
   EMPhotonStandard(const G4String& name = "Photon-Standard");
   virtual ~EMPhotonStandard();
  
 protected:
   void ConstructParticle() {};
   void ConstructProcess();
};

#endif // EMPHOTONSTANDARD_HH








