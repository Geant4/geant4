// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EnergyLossMessenger.hh,v 1.3 2000-11-09 15:52:23 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Class Description:
//  This is a messenger class to interface to exchange information
//  between G4VEnergyLoss and UI.
//        
//  /process/eLoss/   directory
//
//   Commands :
//
//    rndmStep *     Randomize the proposed step by eLoss (false/true)
//    fluct *        Switch true/false the energy loss fluctuations
//    subsec *       Switch true/false the subcutoff generation
//    minsubsec *    Set the min. cut for subcutoff delta in range
//    StepFunction * Set the energy loss step limitation parameters
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef G4EnergyLossMessenger_h
#define G4EnergyLossMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4EnergyLossMessenger: public G4UImessenger
{
  public:   // with description
  
    G4EnergyLossMessenger();
   ~G4EnergyLossMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:

    G4UIdirectory*             eLossDirectory;     
    G4UIcmdWithABool*          RndmStepCmd;
    G4UIcmdWithABool*          EnlossFlucCmd;
    G4UIcmdWithABool*          SubSecCmd;
    G4UIcmdWithADoubleAndUnit* MinSubSecCmd;
    G4UIcommand*               StepFuncCmd;
};

#endif

