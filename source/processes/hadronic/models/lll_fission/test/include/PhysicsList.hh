//******************************************************************************
// PhysicsList.hh
//
// This class is a class derived from G4VUserPhysicsList for constructing 
// particles and physical interaction processes.
//
// 1.00 JMV, LLNL, JAN-2007:  First version.
//******************************************************************************
//
#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class PhysicsList: public G4VUserPhysicsList
{
  public:
    PhysicsList();
    ~PhysicsList();

  protected:
    // Construct particle and physics
    void ConstructParticle();
    void ConstructProcess();

    // Define tracking cuts (step length, etc)
    void SetCuts();
   
  protected:
    // these methods Construct particles 
    void ConstructBosons();
    void ConstructLeptons();
    void ConstructMesons();
    void ConstructBaryons();
    void ConstructIons();

  protected:
    // these methods Construct physics processes and register them
    void ConstructGeneral();
    void ConstructInteractions();


};

#endif
