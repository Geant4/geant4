// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst14PhysicsList.hh,v 1.2 1999-06-06 10:57:18 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst14PhysicsList_h
#define Tst14PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class G4LowEnergyPhotoElectric;
class G4LowEnergyCompton;
class G4LowEnergyRayleigh;
class G4LowEnergyGammaConversion;
class G4LowEnergyBremsstrahlung;

class G4eIonisation;
class G4LowEnergyIonisation;
class G4eBremsstrahlung;
class G4eplusAnnihilation;

class Tst14StepCut;

class Tst14DetectorConstruction;
class Tst14PhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst14PhysicsList: public G4VUserPhysicsList
{
  public:
    Tst14PhysicsList( Tst14DetectorConstruction*);
   ~Tst14PhysicsList();

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

    void SetMaxStep(G4double);

  public:   

    G4double MaxChargedStep;

  private:

    G4LowEnergyPhotoElectric* theLowEnergyPhotoElectric;
    G4LowEnergyCompton*   theLowEnergyCompton;
    G4LowEnergyRayleigh*   theLowEnergyRayleigh;
    G4LowEnergyGammaConversion*     theLowEnergyGammaConversion;
    G4LowEnergyBremsstrahlung*     theLowEnergyBremstrahlung;
    
    G4LowEnergyIonisation*         theeminusIonisation;
    G4eBremsstrahlung*     theeminusBremsstrahlung;
    
    G4eIonisation*         theeplusIonisation;
    G4eBremsstrahlung*     theeplusBremsstrahlung;

    Tst14StepCut* theeminusStepCut ;
    Tst14StepCut* theeplusStepCut ;

    G4double cutForGamma;
    G4double cutForElectron;
    G4double cutForProton;

    Tst14DetectorConstruction* pDet;
    Tst14PhysicsListMessenger* physicsListMessenger;
};

#endif













