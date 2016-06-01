// $Id: ExE01PhysicsList.hh,v 1.1 1998/10/14 15:19:38 allison Exp $
// ------------------------------------------------------------
#ifndef ExE01PhysicsList_h
#define ExE01PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class ExE01PhysicsList: public G4VUserPhysicsList
{
  public:
    ExE01PhysicsList();
    virtual ~ExE01PhysicsList();

  protected:
    // Construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();

    // 
    virtual void SetCuts(G4double cut);
    
  protected:
    // these methods Construct particles 
    virtual void ConstructBosons();
    virtual void ConstructLeptons();

  protected:
  // these methods Construct physics processes and register them
    virtual void ConstructEM();

};

#endif



