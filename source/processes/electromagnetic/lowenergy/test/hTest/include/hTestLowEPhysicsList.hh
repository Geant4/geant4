#ifndef hTestLowEPhysicsList_h
#define hTestLowEPhysicsList_h 1

//---------------------------------------------------------------------------
//
// ClassName:   hTestLowEPhysicsList
//  
// Description: LowEnergy EM processes list
//
// Authors:     V.Ivanchenko 29/03/01
//
// Modified:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestVEMPhysicsList.hh"

// G4 inludes

class hTestLowEPhysicsList:  public hTestVEMPhysicsList
{
public:
  hTestLowEPhysicsList();
  ~hTestLowEPhysicsList() {};

protected:
  void ConstructProcess();
  
private:

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif


