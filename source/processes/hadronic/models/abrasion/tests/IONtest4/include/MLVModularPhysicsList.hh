#ifndef MLVModularPhysicsList_h
#define MLVModularPhysicsList_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "G4VModularPhysicsList.hh"
#include "globals.hh"
////////////////////////////////////////////////////////////////////////////////
//
class MLVModularPhysicsList: public G4VModularPhysicsList
{

public:

  void CleanAllPhysics () {
    G4PhysConstVector::iterator itr;
    for (itr = physicsVector->begin(); itr!= physicsVector->end(); ++itr) {
      delete (*itr);
    }
    physicsVector->clear();
  };
  
};
////////////////////////////////////////////////////////////////////////////////
#endif
