// $Id: Tst05PhysicsList.hh,v 1.4 2000-02-25 16:56:41 gcosmo Exp $
// ------------------------------------------------------------

#ifndef Tst05PhysicsList_h
#define Tst05PhysicsList_h 1

#include "globals.hh"
#include "G4VUserPhysicsList.hh"

class Tst05PhysicsList: public G4VUserPhysicsList
{
  public:
    Tst05PhysicsList();
    virtual ~Tst05PhysicsList();

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



