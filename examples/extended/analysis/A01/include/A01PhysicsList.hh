// $Id: A01PhysicsList.hh,v 1.1 2002-11-13 07:18:49 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
#ifndef A01PhysicsList_h
#define A01PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class A01PhysicsList: public G4VModularPhysicsList
{
public:
  A01PhysicsList();
  virtual ~A01PhysicsList();

public:
  // SetCuts()
  virtual void SetCuts();


};


#endif



