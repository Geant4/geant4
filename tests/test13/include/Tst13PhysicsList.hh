// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst13PhysicsList.hh,v 1.4 1999-12-15 14:54:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef Tst13PhysicsList_h
#define Tst13PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class Tst13PhysicsList: public G4VUserPhysicsList
{
  public:
    Tst13PhysicsList();
    virtual ~Tst13PhysicsList();

  protected:
    // Construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();

    // 
    virtual void SetCuts();
    
  protected:
  // these methods Construct physics processes and register them
    virtual void ConstructGeneral();
    virtual void ConstructEM();
    virtual void ConstructHad();
    virtual void ConstructLeptHad();
 //
    void  ConstructAllBosons();
    void  ConstructAllLeptons();
    void  ConstructAllMesons();
    void  ConstructAllBaryons();
    void  ConstructAllIons();
    void  ConstructAllShortLiveds();

};

#endif



