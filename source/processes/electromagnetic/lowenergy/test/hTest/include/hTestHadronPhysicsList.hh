#ifndef hTestHadronPhysicsList_h
#define hTestHadronPhysicsList_h 1

//---------------------------------------------------------------------------
//
// ClassName:   hTestVHadronPhysicsList
//  
// Description: Standard HARP build Hadron Physics List for Geant4
//
// Authors:     V.Ivanchenko 29/03/01
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "hTestVHadronPhysicsList.hh"

class hTestHadronPhysicsList:  public hTestVHadronPhysicsList
{
public:
  hTestHadronPhysicsList() {};
  ~hTestHadronPhysicsList() {};

protected:
  void ConstructProcess();
  
private:


};

#endif


