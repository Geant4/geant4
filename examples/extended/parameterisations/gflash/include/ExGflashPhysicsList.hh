#ifndef ExGflashPhysicsList_h
#define ExGflashPhysicsList_h 1
using namespace std;

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class ExGflashPhysicsList: public G4VUserPhysicsList
{
public:
  ExGflashPhysicsList();
  virtual ~ExGflashPhysicsList();
  
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
  virtual void ConstructIons();
  
protected:
  // these methods Construct physics processes and register them
  void AddParameterisation();

  virtual void ConstructGeneral();
  virtual void ConstructEM();

  virtual void AddTransportation();
};

#endif



