#ifndef hTestStanPhysicsList_h
#define hTestStanPhysicsList_h 1

//---------------------------------------------------------------------------
//
// ClassName:   hTestStanPhysicsList
//  
// Description: Standard EM physics list
//
// Authors:     V.Ivanchenko 29/03/01
//
// Modified:
//
//----------------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestVEMPhysicsList.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestStanPhysicsList:  public hTestVEMPhysicsList
{
public:
  hTestStanPhysicsList();
  ~hTestStanPhysicsList() {};


protected:
  void ConstructProcess();
  
private:

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif


