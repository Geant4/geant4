//-----------------------------------------------------------
// Renamed of the ufficial Physics List.
//-----------------------------------------------------------
#ifndef CCalPhysicsList_h
#define CCalPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

#define kNuCut  5*m

class CCalPhysicsList: public G4VModularPhysicsList
{
public:
  CCalPhysicsList();
  virtual ~CCalPhysicsList();
  
public:
  // SetCuts() 
  virtual void SetCuts();

};

// 2002 by J.P. Wellisch

#endif



