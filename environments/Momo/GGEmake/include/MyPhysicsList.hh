// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyPhysicsList.hh,v 1.2 1998/11/12 10:50:06 yhajime Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
//
// ExN01PhysicsList
//  Construct/define particles and physics processes
//
//  Particle defined in ExampleN01
//    geantino
//
//  Process defined in ExampleN01
//    transportation
//

#ifndef MyPhysicsList_h
#define MyPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class MyPhysicsList: public G4VUserPhysicsList
{
  public:
    MyPhysicsList();
    ~MyPhysicsList();

  protected:
    // Construct particle and physics process
    void ConstructParticle();
    void ConstructProcess();
    void SetCuts(G4double cut);

};

#endif




