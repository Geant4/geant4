#ifndef Tst01PhysicsList_h
#define Tst01PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class Tst01PhysicsList: public G4VModularPhysicsList
{
public:
  Tst01PhysicsList();
  virtual ~Tst01PhysicsList();
  
public:
  // SetCuts() 
  virtual void SetCuts();

};


#endif



