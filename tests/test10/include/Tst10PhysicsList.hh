// $Id: Tst10PhysicsList.hh,v 1.1 1999-01-08 16:35:32 gunter Exp $
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
    virtual void SetCuts(G4double aCut);
    
  protected:
    // these methods Construct particles 
    virtual void ConstructBosons();

  protected:
  // these methods Construct physics processes and register them
    virtual void ConstructGeneral();
    virtual void ConstructEM();

};

#endif



