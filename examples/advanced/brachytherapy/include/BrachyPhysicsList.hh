//    **********************************
//    *                                *
//    *      BrachyPhysicsList.hh      *
//    *                                *
//    **********************************

#ifndef BrachyPhysicsList_h
#define BrachyPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4LowEnergyIonisation;
class G4LowEnergyPhotoElectric;
class G4LowEnergyBremsstrahlung;

class BrachyPhysicsList: public G4VUserPhysicsList
{
  public:
    BrachyPhysicsList();
   ~BrachyPhysicsList();

  protected:
    // Construct particle and physics
    void ConstructParticle();
    void ConstructProcess();
 
    void SetCuts();
  
  public: 
    // Set Cuts
    void SetGammaCut(G4double);
    void SetElectronCut(G4double);
    void SetPositronCut(G4double);
    
    void SetGammaLowLimit(G4double);
    void SetElectronLowLimit(G4double);
    void SetGELowLimit(G4double);
    void SetLowEnSecPhotCut(G4double);
    void SetLowEnSecElecCut(G4double);
    
  private:
    
    G4double cutForGamma;
    G4double cutForElectron;
    G4double cutForPositron;
    
  protected:
    // these methods Construct particles 
    void ConstructBosons();
    void ConstructLeptons();
    
  protected:
  // these methods Construct physics processes and register them
    void ConstructGeneral();
    void ConstructEM();

  private:
  G4LowEnergyIonisation*  loweIon;
  G4LowEnergyPhotoElectric* lowePhot;
  G4LowEnergyBremsstrahlung* loweBrem;
  
};

#endif



