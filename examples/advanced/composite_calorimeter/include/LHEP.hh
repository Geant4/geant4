#ifndef LHEP_h
#define LHEP_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

#define kNuCut  5*m

class LHEP: public G4VModularPhysicsList
{
public:
  LHEP();
  virtual ~LHEP();
  
public:
  // SetCuts() 
  virtual void SetCuts();

};

// 2002 by J.P. Wellisch

#endif



