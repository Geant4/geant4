#ifndef hTestHadronPhysicsList1_h
#define hTestHadronPhysicsList1_h 1

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
// ClassName:   hTestVHadronPhysicsList1
//  
// Description: hTest Hadron Physics List for Geant4 without ions 
//              and without short lived fragments
//
// Authors:   07.04.01  V.Ivanchenko 
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "hTestVHadronPhysicsList.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestHadronPhysicsList1:  public hTestVHadronPhysicsList
{
public:
  hTestHadronPhysicsList1() {};
  ~hTestHadronPhysicsList1() {};

protected:
  void ConstructProcess();
  
private:

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif


