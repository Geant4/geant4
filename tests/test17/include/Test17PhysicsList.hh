// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Class Description:
// The list of particles and processes are defined in this class.
// Class Description - end
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Test17PhysicsList_h
#define Test17PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class Test17DetectorConstruction;
class Test17PhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Test17PhysicsList: public G4VUserPhysicsList
{
public: // Without description

    Test17PhysicsList( Test17DetectorConstruction*);
   ~Test17PhysicsList();

  protected:
    // Construct particle and physics
    void ConstructParticle();
    void ConstructProcess();
 
    void SetCuts();

  protected:
    // these methods Construct particles 
    void ConstructBosons();
    void ConstructLeptons();
    void ConstructMesons();
    void ConstructBarions();
    void ConstructIons();

  protected:
  // these methods Construct physics processes and register them
    void ConstructGeneral();
    void ConstructEM();
    
  public: // Without description

    void SetGammaCut(G4double);
    void SetElectronCut(G4double);
    void SetProtonCut(G4double);
    void SetCutsByEnergy(G4double);
    void GetRange(G4double);
    void SetMaxStep(G4double);

  private:

    G4double cutForGamma;
    G4double cutForElectron;
    G4double cutForProton;
    
    G4double MaxChargedStep;    

    Test17DetectorConstruction* pDet;
    Test17PhysicsListMessenger* physicsListMessenger;
};

#endif



