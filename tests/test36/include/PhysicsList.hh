#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class StepMax;
class PhysicsListMessenger;
class G4VPhysicsConstructor;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysicsList: public G4VModularPhysicsList
{
  public:
    PhysicsList();
   ~PhysicsList();

    void ConstructParticle();
    void ConstructProcess();
    
    void AddPhysicsList(const G4String& name);
    void AddStepMax();
    
    void SetCuts();
    void SetCutForGamma(G4double);
    void SetCutForElectron(G4double);
    void SetCutForPositron(G4double);

  private:
    G4double cutForGamma;
    G4double cutForElectron;
    G4double cutForPositron;
    G4double currentDefaultCut;
    
    G4VPhysicsConstructor*  emPhysicsList;
    G4String emName;
    
    StepMax* stepMaxProcess;
    
    PhysicsListMessenger* pMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

