// $Id: TstDrawVox01PhysicsList.hh,v 1.1 1999-07-28 17:56:44 graignac Exp $
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This class is a class derived from G4VUserPhysicsList
//      for constructing all particles and processes.
//
//	History
//        first version              10  Jan. 1998 by H.Kurashige
// ------------------------------------------------------------
#ifndef TstDrawVox01PhysicsList_h
#define TstDrawVox01PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class TstDrawVox01PhysicsList: public G4VUserPhysicsList
{
  public:
    TstDrawVox01PhysicsList();
    virtual ~TstDrawVox01PhysicsList();

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



