// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst17PhysicsList.hh,v 1.1 1999-11-30 18:01:51 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst17PhysicsList_h
#define Tst17PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class G4MultipleScattering;

class G4eIonisation;
class G4eBremsstrahlung;
class G4eplusAnnihilation;

class G4LowEnergyPhotoElectric;
class G4LowEnergyCompton;
class G4LowEnergyGammaConversion;
class G4LowEnergyRayleigh;

class G4LowEnergyIonisation;
class G4LowEnergyBremsstrahlung;
    
class Tst17DetectorConstruction;
class Tst17PhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst17PhysicsList: public G4VUserPhysicsList
{
  public:
    Tst17PhysicsList(Tst17DetectorConstruction*);
   ~Tst17PhysicsList();

  protected:
    // Construct particle and physics
    void ConstructParticle();
    void ConstructProcess();
 
    void SetCuts();

  protected:
    // these methods Construct particles 
    void ConstructBosons();
    void ConstructLeptons();

  protected:
  // these methods Construct physics processes and register them
    void ConstructGeneral();
    void ConstructEM();
    
public:
  
  void SetGammaCut(G4double);
  void SetElectronCut(G4double);
  void SetProtonCut(G4double);
  void SetCutsByEnergy(G4double);
  void GetRange(G4double);

  void SetGammaLowLimit(G4double);
  void SetElectronLowLimit(G4double);
  void SetGELowLimit(G4double);
  void SetLowEnSecPhotCut(G4double);
  void SetLowEnSecElecCut(G4double);

public:   
  
  G4double MaxChargedStep;
  
private:
  
  G4LowEnergyPhotoElectric*    theLEPhotoElectric;
  G4LowEnergyCompton*          theLECompton;
  G4LowEnergyGammaConversion*  theLEGammaConversion;
  G4LowEnergyRayleigh*         theLERayleigh;
  
  G4MultipleScattering*  theeminusMultipleScattering;
  G4LowEnergyIonisation*       theLEIonisation;
  G4LowEnergyBremsstrahlung*   theLEBremsstrahlung;
  
  G4MultipleScattering*  theeplusMultipleScattering;
  G4eIonisation*         theeplusIonisation;
  G4eBremsstrahlung*     theeplusBremsstrahlung;
  G4eplusAnnihilation*   theeplusAnnihilation;
  
  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForProton;
  
  Tst17DetectorConstruction* pDet;
  Tst17PhysicsListMessenger* physicsListMessenger;
};

#endif
















