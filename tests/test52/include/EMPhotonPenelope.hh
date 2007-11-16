
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
//  The considered processes are G4PenelopePhotoElectric, G4PenelopeCompton
//  G4PenelopeGammaConversion and G4PenelopeRayleigh 
//

#ifndef EMPHOTONPENELOPE_HH
#define EMPHOTONPENELOPE_HH 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"


class EMPhotonPenelope : public G4VPhysicsConstructor {

 public: 
   EMPhotonPenelope(const G4String& name = "Photon-Penelope"); 
   virtual ~EMPhotonPenelope();
  
 protected:
   void ConstructParticle() {};
   void ConstructProcess();
};

#endif // EMPHOTONPENELOPE_HH

