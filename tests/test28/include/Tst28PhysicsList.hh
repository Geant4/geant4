#ifndef Tst28PhysicsList_h
#define Tst28PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class Tst28PhysicsList: public G4VUserPhysicsList
{
  public:
    Tst28PhysicsList();
    virtual ~Tst28PhysicsList();

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
    virtual void ConstructLeptHad();
 //
    void  ConstructAllBosons();
    void  ConstructAllLeptons();
    void  ConstructAllMesons();
    void  ConstructAllBaryons();
    void  ConstructAllIons();
    void  ConstructAllShortLiveds();

};

#endif



