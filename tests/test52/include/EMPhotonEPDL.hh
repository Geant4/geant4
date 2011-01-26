
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
//  The considered processes are G4LivermorePhotoElectric, G4LivermoreCompton
//  G4LivermoreGammaConversion and G4LivermoreRayleigh (Low Energy Package, 
//  EEDL libraries)
//

#ifndef EMPHOTONEPDL_HH
#define EMPHOTONEPDL_HH 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"


class EMPhotonEPDL : public G4VPhysicsConstructor {

 public: 
   EMPhotonEPDL(const G4String& name = "Photon-EPDL"); 
   virtual ~EMPhotonEPDL();
  
 protected:
   void ConstructParticle() {};
   void ConstructProcess();
};

#endif // EMPHOTONEPDL_HH








