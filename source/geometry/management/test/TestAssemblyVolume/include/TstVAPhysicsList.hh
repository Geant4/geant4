// $Id: TstVAPhysicsList.hh,v 1.2 2001-02-01 21:25:36 radoone Exp $
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This class is a class derived from G4VUserPhysicsList
//      for constructing all particles and processes.
//
//	History
//        first version              10  Jan. 1998 by H.Kurashige
// ------------------------------------------------------------
#ifndef TstVAPhysicsList_h
#define TstVAPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class TstVAPhysicsList: public G4VUserPhysicsList
{
  public:
    TstVAPhysicsList();
    virtual ~TstVAPhysicsList();

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



