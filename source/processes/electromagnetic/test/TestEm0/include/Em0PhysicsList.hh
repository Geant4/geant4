// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em0PhysicsList.hh,v 1.3 1999-05-10 16:44:38 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em0PhysicsList_h
#define Em0PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class Em0DetectorConstruction;
class Em0PhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em0PhysicsList: public G4VUserPhysicsList
{
  public:
    Em0PhysicsList(Em0DetectorConstruction*);
   ~Em0PhysicsList();

  private:
    Em0DetectorConstruction* pDet;
    Em0PhysicsListMessenger* pMes;

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
    void SetCutsByEnergy(G4double);
    void GetRange(G4double);
    void SetECut(G4double);
  
  private:
    G4double cutForGamma;
    G4double cutForE;
};

#endif



