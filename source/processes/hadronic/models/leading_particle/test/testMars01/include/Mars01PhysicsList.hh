// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Mars01PhysicsList.hh,v 1.1 2001-12-13 14:58:43 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Mars01PhysicsList_h
#define Mars01PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Mars01PhysicsList: public G4VUserPhysicsList
{
public:
  Mars01PhysicsList();
  ~Mars01PhysicsList();
  
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
  void ConstructBaryons();
  
protected:
  // these methods Construct physics processes and register them
  void ConstructGeneral();
  void ConstructEM();
  void ConstructMarsProc();
};

#endif

 
