#ifndef hTestHadronPhysicsList1_h
#define hTestHadronPhysicsList1_h 1

//---------------------------------------------------------------------------
//
// ClassName:   hTestVHadronPhysicsList1
//  
// Description: HARP build Hadron Physics List for Geant4 without ions 
//              and without short lived fragments
//
// Authors:     V.Ivanchenko 29/03/01
//
// Modified:
//
//----------------------------------------------------------------------------
//

// HARP inludes
#include "hTestVHadronPhysicsList.h"

// G4 inludes

class hTestHadronPhysicsList1:  public hTestVHadronPhysicsList
{
public:
  hTestHadronPhysicsList1() {};
  ~hTestHadronPhysicsList1() {};

protected:
  void ConstructProcess();
  
private:


};

#endif


