#ifndef Tst02PhysicsList_h
#define Tst02PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class Tst02PhysicsList: public G4VModularPhysicsList
{
public:
  Tst02PhysicsList();
  virtual ~Tst02PhysicsList();
  
public:
  // SetCuts() 
  virtual void SetCuts();


};


#endif



