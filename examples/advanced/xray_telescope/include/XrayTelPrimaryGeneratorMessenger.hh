// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelPrimaryGeneratorMessenger.hh             *
// * -------                                                            *
// *                                                                    *
// * Version:           0.4                                             *
// * Date:              06/11/00                                        *
// * Author:            R Nartallo                                      *
// * Organisation:      ESA/ESTEC, Noordwijk, THe Netherlands           *
// *                                                                    *
// **********************************************************************
// 
// CHANGE HISTORY
// --------------
//
// 06.11.2000 R.Nartallo
// - First implementation of X-ray Telescope advanced example.
// - Based on Chandra and XMM models
//
//
// **********************************************************************

#ifndef XrayTelPrimaryGeneratorMessenger_h
#define XrayTelPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class XrayTelPrimaryGeneratorAction;
class G4UIcmdWithAString;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayTelPrimaryGeneratorMessenger: public G4UImessenger
{
public:
  XrayTelPrimaryGeneratorMessenger(XrayTelPrimaryGeneratorAction*);
  ~XrayTelPrimaryGeneratorMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  XrayTelPrimaryGeneratorAction*    XrayTelAction; 
  G4UIcmdWithAString*          RndmCmd;
  G4UIcmdWithAString*          ErndmCmd;
  G4UIcmdWithADoubleAndUnit*   SetRmin;
  G4UIcmdWithADoubleAndUnit*   SetRmax;
  G4UIcmdWithADoubleAndUnit*   SetTmax;
};

#endif

