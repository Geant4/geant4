
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef exrdmPhysListGeneral_h
#define exrdmPhysListGeneral_h 1

#include "G4VPhysicsConstructor.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class exrdmPhysListGeneral : public G4VPhysicsConstructor
{
public: 
  exrdmPhysListGeneral(const G4String& name = "general");
  virtual ~exrdmPhysListGeneral();

public: 
  // This method is dummy for physics
  virtual void ConstructParticle() {};
 
  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type 
  virtual void ConstructProcess();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif








