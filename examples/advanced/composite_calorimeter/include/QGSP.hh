#ifndef TQGSP_h
#define TQGSP_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

template<class T>
class TQGSP: public T
{
public:
  TQGSP();
  virtual ~TQGSP();
  
public:
  // SetCuts() 
  virtual void SetCuts();

};
#include "QGSP.icc"
typedef TQGSP<G4VModularPhysicsList> QGSP;

// 2002 by J.P. Wellisch

#endif



