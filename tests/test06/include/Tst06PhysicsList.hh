// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst06PhysicsList.hh,v 1.1 1999-01-08 16:35:04 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst06PhysicsList_h
#define Tst06PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst06PhysicsList: public G4VUserPhysicsList
{
  public:
    Tst06PhysicsList();
   ~Tst06PhysicsList();

  protected:
    // Construct particles
    void ConstructParticle();
    void ConstructBosons();
    void ConstructLeptons();
    void ConstructMesons();
    void ConstructBarions(); 
     
  public:
    void SetCuts(G4double);
    void SetGammaCut(G4double);
    void SetElectronCut(G4double);
    void SetProtonCut(G4double);           
        
  protected:
    // Construct processes and register them
    void ConstructProcess();  
    void ConstructGeneral();
    void ConstructEM();
    
  private:

    G4double cutForGamma;
    G4double cutForElectron;
    G4double cutForProton;

};

#endif

