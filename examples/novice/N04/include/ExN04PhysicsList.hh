#ifndef ExN04PhysicsList_h
#define ExN04PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class ExN04PhysicsList: public G4VUserPhysicsList
{
public:
  ExN04PhysicsList();
  virtual ~ExN04PhysicsList();
  
protected:
  // Construct particle and physics
  virtual void ConstructParticle();
  virtual void ConstructProcess();
  
  // 
  virtual void SetCuts();
  
protected:
  // these methods Construct physics processes and register them
  virtual void ConstructGeneral();
  virtual void ConstructEM();
  virtual void ConstructHad();

  // these methods Construct all particles in each category
  virtual void ConstructAllBosons();
  virtual void ConstructAllLeptons();
  virtual void ConstructAllMesons();
  virtual void ConstructAllBarions();
  virtual void ConstructAllIons();
  virtual void ConstructAllShortLiveds();
};

#endif



