//  ExamplePhysicsList class description
// ------------------------------------------------------------
#ifndef ExamplePhysicsList_h
#define ExamplePhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class ExamplePhysicsList: public G4VUserPhysicsList
{
  public:
    ExamplePhysicsList();
    virtual ~ExamplePhysicsList();

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
    virtual void ConstructGeneral();
    virtual void ConstructEM();
    virtual void ConstructHad();
    virtual void ConstructLeptHad();

};

#endif







