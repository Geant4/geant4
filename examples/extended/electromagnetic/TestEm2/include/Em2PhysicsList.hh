// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em2PhysicsList.hh,v 1.3 2000-04-17 11:25:52 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em2PhysicsList_h
#define Em2PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class Em2PhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em2PhysicsList: public G4VUserPhysicsList
{
  public:
    Em2PhysicsList();
   ~Em2PhysicsList();

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
    
    Em2PhysicsListMessenger* pMessenger;         
};

#endif

