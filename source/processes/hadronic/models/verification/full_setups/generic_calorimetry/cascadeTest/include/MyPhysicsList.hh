#ifndef MyPhysicsList_h
#define MyPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

#define kNuCut  5*m

class MyPhysicsList: public G4VModularPhysicsList {

public:

  MyPhysicsList();
  virtual ~MyPhysicsList();
  
  // SetCuts() 
  virtual void SetCuts();

};

#endif



