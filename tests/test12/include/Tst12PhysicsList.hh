// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst12PhysicsList.hh,v 1.3 1999-10-03 10:07:23 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef Tst12PhysicsList_h
#define Tst12PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class Tst12PhysicsList: public G4VUserPhysicsList
{
  public:
    Tst12PhysicsList();
    virtual ~Tst12PhysicsList();

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



