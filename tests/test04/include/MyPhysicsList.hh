// $Id: MyPhysicsList.hh,v 1.1 1999-01-08 16:34:46 gunter Exp $
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This class is a class derived from G4VUserPhysicsList
//      for constructing all particles and processes.
//
//	History
//        first version              10  Jan. 1998 by H.Kurashige
// ------------------------------------------------------------
#ifndef MyPhysicsList_h
#define MyPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class MyPhysicsList: public G4VUserPhysicsList
{
  public:
    MyPhysicsList();
    virtual ~MyPhysicsList();

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



