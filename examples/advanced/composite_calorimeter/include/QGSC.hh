#ifndef TQGSC_h
#define TQGSC_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

template<class T>
class TQGSC: public T
{
public:
  TQGSC();
  virtual ~TQGSC();
  
public:
  // SetCuts() 
  virtual void SetCuts();

};
#include "QGSC.icc"
typedef TQGSC<G4VModularPhysicsList> QGSC;

// 2002 by J.P. Wellisch

#endif



