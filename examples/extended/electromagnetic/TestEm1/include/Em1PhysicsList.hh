// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em1PhysicsList.hh,v 1.1 1999-10-11 13:07:37 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em1PhysicsList_h
#define Em1PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class Em1DetectorConstruction;
class Em1PhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em1PhysicsList: public G4VUserPhysicsList
{
  public:
    Em1PhysicsList(Em1DetectorConstruction*);
   ~Em1PhysicsList();

  protected:
    // Construct particles
    void ConstructParticle();
    void ConstructBosons();
    void ConstructLeptons();
    void ConstructMesons();
    void ConstructBarions(); 
     
  public:
    void SetCuts();
    void SetGammaCut(G4double);
    void SetElectronCut(G4double);
    void SetProtonCut(G4double);           
    void SetCutsByEnergy(G4double);
    void GetRange(G4double);  
        
  protected:
    // Construct processes and register them
    void ConstructProcess();  
    void ConstructGeneral();
    void ConstructEM();
    
  private:
    Em1DetectorConstruction* pDet;
    Em1PhysicsListMessenger* pMes;    

    G4double cutForGamma;
    G4double cutForElectron;
    G4double cutForProton;

};

#endif

