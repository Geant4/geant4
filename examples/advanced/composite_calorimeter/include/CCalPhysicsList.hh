///////////////////////////////////////////////////////////////////////////////
// File: CCalPhysicsList.hh
// Description: CCalPhysicsList provides Physics list for the simulation
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalPhysicsList_h
#define CCalPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class CCalPhysicsList: public G4VModularPhysicsList {

public:
  CCalPhysicsList();
  virtual ~CCalPhysicsList();
  
protected:
  //  SetCuts()
  virtual void SetCuts();
  
};

#endif



