// $Id: ExE02PhysicsList.hh,v 1.2 1999/04/17 04:50:02 kurasige Exp $
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This class is a class derived from G4VUserPhysicsList
//      for constructing all particles and processes.
//
//	History
//        first version              10  Jan. 1998 by H.Kurashige
// ------------------------------------------------------------
#ifndef ExE02PhysicsList_h
#define ExE02PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class ExE02PhysicsList: public G4VUserPhysicsList
{
  public:
    ExE02PhysicsList();
    virtual ~ExE02PhysicsList();

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



