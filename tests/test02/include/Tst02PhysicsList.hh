#ifndef Tst02PhysicsList_h
#define Tst02PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class Tst02PhysicsList: public G4VUserPhysicsList
{
public:
  Tst02PhysicsList();
  virtual ~Tst02PhysicsList();
  
protected:
  // Construct particle and physics
  virtual void ConstructParticle();
  virtual void ConstructProcess();
  
  // 
  virtual void SetCuts();
  
protected:
  // these methods Construct physics processes and register them
  void AddParameterisation();
  virtual void ConstructGeneral();
  virtual void ConstructEM();
  virtual void ConstructHad();
  virtual void ConstructLeptHad();

  // these methods Construct all particles in each category
  virtual void ConstructAllBosons();
  virtual void ConstructAllLeptons();
  virtual void ConstructAllMesons();
  virtual void ConstructAllBarions();
  virtual void ConstructAllIons();
  virtual void ConstructAllShortLiveds();
};

#endif



