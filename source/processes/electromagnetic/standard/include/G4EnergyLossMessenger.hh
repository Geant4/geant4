// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EnergyLossMessenger.hh,v 1.3 2000-05-23 14:39:27 urban Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef G4EnergyLossMessenger_h
#define G4EnergyLossMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIcommand;
class G4UIcmdWithABool;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4EnergyLossMessenger: public G4UImessenger
{
  public:
  
    G4EnergyLossMessenger();
   ~G4EnergyLossMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
      
    G4UIcmdWithABool*          RndmStepCmd;
    G4UIcmdWithABool*          EnlossFlucCmd;
    G4UIcmdWithABool*          SubSecCmd;
    G4UIcommand*               StepFuncCmd;
};

#endif

