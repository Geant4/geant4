
#ifndef LXePhysicsList_h
#define LXePhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class LXePhysicsList: public G4VModularPhysicsList
{
public:
  LXePhysicsList();
  virtual ~LXePhysicsList();
  
public:
  // SetCuts() 
  virtual void SetCuts();


};


#endif



