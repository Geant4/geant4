#ifndef ExN04PhysicsList_h
#define ExN04PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class ExN04PhysicsList: public G4VModularPhysicsList
{
public:
  ExN04PhysicsList();
  virtual ~ExN04PhysicsList();
  
public:
  // SetCuts() 
  virtual void SetCuts();


};


#endif



