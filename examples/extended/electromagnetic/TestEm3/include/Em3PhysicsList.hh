// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em3PhysicsList.hh,v 1.3 2000-04-17 12:06:23 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em3PhysicsList_h
#define Em3PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class Em3PhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em3PhysicsList: public G4VUserPhysicsList
{
  public:
    Em3PhysicsList();
   ~Em3PhysicsList();

  protected:
    // Construct particle
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
  // Construct physics processes and register them
    void ConstructProcess();  
    void ConstructGeneral();
    void ConstructEM();
    
  private:
    G4double cutForGamma;
    G4double cutForElectron; 
    G4double cutForProton;
    G4double currentDefaultCut;
    
    Em3PhysicsListMessenger* pMessenger;             
};

#endif



