
#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class StepMax;
class PhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysicsList: public G4VModularPhysicsList
{
  public:
    PhysicsList();
   ~PhysicsList();

    void ConstructParticle();
    
    void SetCuts();
    void SetCutForGamma(G4double);
    void SetCutForElectron(G4double);
    void SetCutForPositron(G4double);        
        
    void AddPhysicsList(const G4String& name);
    void ConstructProcess();
    
    void AddDecay();
    void AddStepMax();       

  private:
    G4double cutForGamma;
    G4double cutForElectron;
    G4double cutForPositron;
       
    G4String                             emName;
    G4VPhysicsConstructor*               emPhysicsList;    
    
    PhysicsListMessenger* pMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

