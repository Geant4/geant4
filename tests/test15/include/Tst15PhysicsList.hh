#ifndef Tst15PhysicsList_h
#define Tst15PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class Tst15PhysicsList: public G4VUserPhysicsList
{
public:
  Tst15PhysicsList();
  virtual ~Tst15PhysicsList();
  
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
  virtual void ConstructAllBaryons();
  virtual void ConstructAllIons();
  virtual void ConstructAllShortLiveds();
};

#endif



