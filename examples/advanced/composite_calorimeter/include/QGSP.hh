#ifndef QGSP_h
#define QGSP_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

#define kNuCut  5*m

class QGSP: public G4VModularPhysicsList
{
public:
  QGSP();
  virtual ~QGSP();
  
public:
  // SetCuts() 
  virtual void SetCuts();

};

// 2002 by J.P. Wellisch

#endif



