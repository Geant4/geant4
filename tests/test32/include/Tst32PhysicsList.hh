// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst32PhysicsList.hh,v 1.1 2002-06-13 12:16:35 jwellisc Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst32PhysicsList_h
#define Tst32PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst32PhysicsList: public G4VUserPhysicsList
{
public:
  Tst32PhysicsList();
  ~Tst32PhysicsList();
  
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

 
