#ifndef TLHEP_h
#define TLHEP_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

template<class T>
class TLHEP: public T
{
public:
  TLHEP();
  virtual ~TLHEP();
  
public:
  // SetCuts() 
  virtual void SetCuts();

};
#include "LHEP.icc"
typedef TLHEP<G4VModularPhysicsList> LHEP;

// 2002 by J.P. Wellisch

#endif



