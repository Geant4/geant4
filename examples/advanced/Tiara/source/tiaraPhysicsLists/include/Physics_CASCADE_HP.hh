#ifndef Physics_CASCADE_HP_h
#define Physics_CASCADE_HP_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

#define kNuCut  5*m

class Physics_CASCADE_HP: public G4VModularPhysicsList
{
public:
  Physics_CASCADE_HP();
  virtual ~Physics_CASCADE_HP();
  
public:
  // SetCuts() 
  virtual void SetCuts();

};

// 2002 by J.P. Wellisch

#endif



