// $Id: Tst05PhysicsList.hh,v 1.1 1999-01-08 16:34:50 gunter Exp $
// ------------------------------------------------------------
#ifndef Tst05PhysicsList_h
#define Tst05PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

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
    virtual void SetCuts(G4double aCut);
    
  protected:
    // these methods Construct particles 
    virtual void ConstructBosons();
    virtual void ConstructLeptons();
    virtual void ConstructMesons();
    virtual void ConstructBarions();
    virtual void ConstructIons();

  protected:
  // these methods Construct physics processes and register them
    virtual void ConstructGeneral();
    virtual void ConstructEM();
    virtual void ConstructHad();
    virtual void ConstructLeptHad();

};

#endif



