#ifndef LHEP_BIC_HP_h
#define LHEP_BIC_HP_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

#define kNuCut  5*m

class LHEP_BIC_HP: public G4VModularPhysicsList
{
public:
  LHEP_BIC_HP();
  virtual ~LHEP_BIC_HP();
  
public:
  // SetCuts() 
  virtual void SetCuts();

};

// 2002 by J.P. Wellisch

#endif



