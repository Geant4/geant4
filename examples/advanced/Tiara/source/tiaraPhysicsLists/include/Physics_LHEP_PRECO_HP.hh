#ifndef Physics_LHEP_PRECO_HP_h
#define Physics_LHEP_PRECO_HP_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

#define kNuCut  5*m

class Physics_LHEP_PRECO_HP: public G4VModularPhysicsList
{
public:
  Physics_LHEP_PRECO_HP();
  virtual ~Physics_LHEP_PRECO_HP();
  
public:
  // SetCuts() 
  virtual void SetCuts();

};

// 2002 by J.P. Wellisch

#endif



