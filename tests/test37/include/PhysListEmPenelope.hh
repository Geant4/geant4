#ifndef PhysListEmPenelope_h
#define PhysListEmPenelope_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysListEmPenelope : public G4VPhysicsConstructor
{
  public: 
    PhysListEmPenelope(const G4String& name = "Penelope");
   ~PhysListEmPenelope();

  public: 
    // This method is dummy for physics
    void ConstructParticle() {};
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    void ConstructProcess();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif








