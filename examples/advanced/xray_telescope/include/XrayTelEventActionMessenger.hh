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
// * MODULE:            XrayTelEventActionMessenger.hh                  *
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

#ifndef XrayTelEventActionMessenger_h
#define XrayTelEventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class XrayTelEventAction;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayTelEventActionMessenger: public G4UImessenger
{
public:
  XrayTelEventActionMessenger(XrayTelEventAction*);
  ~XrayTelEventActionMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  XrayTelEventAction*   eventAction;   
  G4UIcmdWithAString* DrawCmd;
};

#endif






