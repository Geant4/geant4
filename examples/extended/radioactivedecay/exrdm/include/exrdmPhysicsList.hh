
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef exrdmPhysicsList_h
#define exrdmPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4VPhysicsConstructor;
class exrdmPhysicsListMessenger;
class G4ProductionCuts;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class exrdmPhysicsList: public G4VModularPhysicsList
{
public:
  exrdmPhysicsList();
  ~exrdmPhysicsList();

  void ConstructParticle();

  void SetCuts();
  void SetCutForGamma(G4double);
  void SetCutForElectron(G4double);
  void SetCutForPositron(G4double);

  void SelectPhysicsList(const G4String& name);
  void ConstructProcess();

  void SetTargetCut(G4double val);
  void SetDetectorCut(G4double val);

private:

  // hide assignment operator
  exrdmPhysicsList & operator=(const exrdmPhysicsList &right);
  exrdmPhysicsList(const exrdmPhysicsList&);

  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;

  G4VPhysicsConstructor*  emPhysicsList;
  G4VPhysicsConstructor*  generalPhysicsList;
  G4VPhysicsConstructor*  particleList;
  //  std::vector<G4VPhysicsConstructor*>  hadronPhys;
  G4VPhysicsConstructor*  hadPhysicsList;

  exrdmPhysicsListMessenger* pMessenger;
  G4ProductionCuts* DetectorCuts;
  G4ProductionCuts* TargetCuts;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

