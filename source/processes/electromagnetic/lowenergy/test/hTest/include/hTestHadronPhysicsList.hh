#ifndef hTestHadronPhysicsList_h
#define hTestHadronPhysicsList_h 1

//---------------------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//---------------------------------------------------------------------------
//
// ClassName:   hTestVHadronPhysicsList
//  
// Description: Standard hTest Hadron Physics List for Geant4
//
// Authors:     V.Ivanchenko 29/03/01
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "hTestVHadronPhysicsList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestHadronPhysicsList:  public hTestVHadronPhysicsList
{
public:
  hTestHadronPhysicsList() {};
  ~hTestHadronPhysicsList() {};

protected:
  void ConstructProcess();
  
private:


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif


