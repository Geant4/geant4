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
  virtual void SetCuts();
  
protected:
  // these methods Construct particles 
  virtual void ConstructBosons();
  virtual void ConstructLeptons();
  virtual void ConstructMesons();
  virtual void ConstructBaryons();
  
protected:
  // these methods Construct physics processes and register them
  void AddParameterisation();
  virtual void ConstructGeneral();
  virtual void ConstructEM();

};

#endif



