#ifndef Physics_LHEP_PRECO_h
#define Physics_LHEP_PRECO_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

#define kNuCut  5*m

class Physics_LHEP_PRECO: public G4VModularPhysicsList
{
public:
  Physics_LHEP_PRECO();
  virtual ~Physics_LHEP_PRECO();

  
public:
  // SetCuts() 
  virtual void SetCuts();

};

// 2002 by J.P. Wellisch

#endif



