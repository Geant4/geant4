// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em5PhysicsList.hh,v 1.3 2001-01-09 12:36:37 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em5PhysicsList_h
#define Em5PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class Em5DetectorConstruction;
class Em5PhysicsListMessenger;
class Em5StepCut;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em5PhysicsList: public G4VUserPhysicsList
{
  public:
    Em5PhysicsList( Em5DetectorConstruction*);
   ~Em5PhysicsList();

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

  private:

    G4double cutForGamma;
    G4double cutForElectron;
    G4double cutForProton;
    G4double currentDefaultCut;

    Em5DetectorConstruction* pDet;
    Em5PhysicsListMessenger* physicsListMessenger;
    Em5StepCut* pStepCut;
};

#endif



