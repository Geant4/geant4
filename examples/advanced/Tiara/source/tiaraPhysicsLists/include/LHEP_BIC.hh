#ifndef LHEP_BIC_h
#define LHEP_BIC_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

#define kNuCut  5*m

class LHEP_BIC: public G4VModularPhysicsList
{
public:
  LHEP_BIC();
  virtual ~LHEP_BIC();
  
public:
  // SetCuts() 
  virtual void SetCuts();

};

// 2002 by J.P. Wellisch

#endif



