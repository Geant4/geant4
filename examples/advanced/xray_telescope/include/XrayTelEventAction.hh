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
// * MODULE:            XrayTelEventAction.hh                           *
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

#ifndef XrayTelEventAction_h
#define XrayTelEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class XrayTelEventActionMessenger;
class XrayTelAnalysisManager;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayTelEventAction : public G4UserEventAction
{
public:
  XrayTelEventAction(XrayTelAnalysisManager* = 0);
  ~XrayTelEventAction();

public:
  void BeginOfEventAction(const G4Event* anEvent);
  void EndOfEventAction(const G4Event* anEvent);
    
  void SetDrawFlag(G4String val)  {fDrawFlag = val;};
    
private:
  XrayTelAnalysisManager* fAnalysisManager;
  G4String fDrawFlag;                         // control the drawing of event
  XrayTelEventActionMessenger* fEventMessenger;
};

#endif

    




