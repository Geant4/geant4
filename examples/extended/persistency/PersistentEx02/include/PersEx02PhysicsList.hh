// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: PersEx02PhysicsList.hh,v 1.3 1999/11/29 18:23:32 morita Exp $
// GEANT4 tag $Name: geant4-02-00 $
//

#ifndef PersEx02PhysicsList_h
#define PersEx02PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class PersEx02PhysicsList: public G4VUserPhysicsList
{
  public:
    PersEx02PhysicsList();
    virtual ~PersEx02PhysicsList();

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



