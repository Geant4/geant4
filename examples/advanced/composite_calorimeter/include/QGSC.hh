#ifndef QGSC_h
#define QGSC_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

#define kNuCut  5*m

class QGSC: public G4VModularPhysicsList
{
public:
  QGSC();
  virtual ~QGSC();
  
public:
  // SetCuts() 
  virtual void SetCuts();

};

// 2002 by J.P. Wellisch

#endif



