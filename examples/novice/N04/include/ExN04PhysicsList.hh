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
  virtual void SetCuts(G4double aCut);
  
protected:
  // these methods Construct physics processes and register them
  void AddParameterisation();
  virtual void ConstructGeneral();
  virtual void ConstructEM();
  virtual void ConstructHad();
  virtual void ConstructLeptHad();
};

#endif



