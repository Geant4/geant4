//  Em6PhysicsList.hh

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Em6PhysicsList_h
#define Em6PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"
#include "G4GammaConversionToMuons.hh"

class Em6PhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Em6PhysicsList: public G4VUserPhysicsList
{
  public:
    Em6PhysicsList();
   ~Em6PhysicsList();

  protected:
    // Construct particles
    void ConstructParticle();
    void ConstructBosons();
    void ConstructLeptons();
    void ConstructMesons();
    void ConstructBarions();

  public:
    void SetCuts();
    void SetCutForGamma(G4double);
    void SetCutForElectron(G4double);
    void SetCutForProton(G4double);
    void setGammaToMuPairFac(G4double);

  protected:
    // Construct processes and register them
    void ConstructProcess();
    void ConstructGeneral();
    void ConstructEM();

  private:
    G4double cutForGamma;
    G4double cutForElectron;
    G4double cutForProton;
    G4double currentDefaultCut;
    G4GammaConversionToMuons *theGammaToMuPairProcess;

    Em6PhysicsListMessenger* pMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

