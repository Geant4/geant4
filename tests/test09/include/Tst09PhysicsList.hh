#ifndef Tst09PhysicsList_h
#define Tst09PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class Tst09PhysicsList: public G4VUserPhysicsList
{
public:
  Tst09PhysicsList();
  virtual ~Tst09PhysicsList();
  
protected:
  // Construct particle and physics
  virtual void ConstructParticle();
  virtual void ConstructProcess();
  
  // 
  virtual void SetCuts(G4double aCut);
  
protected:
  // these methods Construct particles 
  virtual void ConstructBosons();
  virtual void ConstructLeptons();
  virtual void ConstructMesons();
  virtual void ConstructBarions();
  
protected:
  // these methods Construct physics processes and register them
  void AddParameterisation();
  virtual void ConstructGeneral();
  virtual void ConstructEM();

  G4double defaultCutValue;
  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForProton;
};

#endif



