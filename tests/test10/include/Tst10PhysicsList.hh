// $Id: Tst10PhysicsList.hh,v 1.2 1999-04-17 08:01:47 kurasige Exp $
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This class is a class derived from G4VUserPhysicsList
//      for constructing all particles and processes.
//
//	History
//        first version              09 Sept. 1998 by S.Magni
// ------------------------------------------------------------
#ifndef Tst10PhysicsList_h
#define Tst10PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class Tst10PhysicsList: public G4VUserPhysicsList
{
  public:
    Tst10PhysicsList();
    virtual ~Tst10PhysicsList();

  protected:
    // Construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();

    // 
    virtual void SetCuts();
    
  protected:
    // these methods Construct particles 
    virtual void ConstructBosons();

  protected:
  // these methods Construct physics processes and register them
    virtual void ConstructGeneral();
    virtual void ConstructEM();

};

#endif



