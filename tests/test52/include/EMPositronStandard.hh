
// ==============================
// PHYSICS PROCESSES:
// ==============================
//  Electromagnetic processes for: 
//      positrons
//     
// 
// ==============================
// COMMENTS:
// ==============================
//  The considered processes are G4MultipleScattering, G4eIonisation,
//  G4eBremsstrahlung and G4eplusAnnihilation
//

#ifndef EMPOSITRONSTANDARD_HH
#define EMPOSITRONSTANDARD_HH 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"


class EMPositronStandard : public G4VPhysicsConstructor {

 public: 
   EMPositronStandard(const G4String& name = "Positron-Standard");
   virtual ~EMPositronStandard();
  
 protected:
   void ConstructParticle() {};
   void ConstructProcess();

   G4double facRange;
};

#endif // EMPOSITRONSTANDARD_HH








