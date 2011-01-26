
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
//  The considered processes are G4MultipleScattering, G4PenelopeIonisation,
//  G4PenelopeBremsstrahlung and G4PenelopeAnnihilation
//

#ifndef EMPOSITRONPENELOPE_HH
#define EMPOSITRONPENELOPE_HH 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"


class EMPositronPenelope : public G4VPhysicsConstructor {

 public: 
   EMPositronPenelope(const G4String& name = "Positron-Penelope");
   virtual ~EMPositronPenelope();
  
 protected:
   void ConstructParticle() {};
   void ConstructProcess();

   G4double facRange;
};

#endif // EMPOSITRONPENELOPE_HH








