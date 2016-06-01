// $Id: ExE03PhysicsList.hh,v 1.1 1998/10/14 15:55:45 allison Exp $
// ------------------------------------------------------------
#ifndef ExE03PhysicsList_h
#define ExE03PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class ExE03PhysicsList: public G4VUserPhysicsList
{
  public:
    ExE03PhysicsList();
    virtual ~ExE03PhysicsList();

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



