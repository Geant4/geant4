// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#ifndef Tst03PhysicsList_h
#define Tst03PhysicsList_h 1

#include "globals.hh"
#include "G4VUserPhysicsList.hh"

class Tst03PhysicsList : public G4VUserPhysicsList
{
  public:
    Tst03PhysicsList();
    virtual ~Tst03PhysicsList();

  protected:
    // Construct particles and processes
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

  protected:
  // these methods Construct physics processes and register them
  virtual void ConstructGeneral();
  virtual void ConstructEM();
  virtual void ConstructOp();

};

#endif /* Tst03PhysicsList_h */




