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
// * MODULE:            XrayTelSteppingAction.cc                        *
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
// - First implementation of xray_telescope SteppingAction
// - Based on Chandra and XMM models
// 
//
// **********************************************************************

#include "G4SteppingManager.hh"

#include "XrayTelSteppingAction.hh"

#ifdef G4ANALYSIS_USE
#include "XrayTelAnalysisManager.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayTelSteppingAction::XrayTelSteppingAction(
 XrayTelAnalysisManager* aAnalysisManager
):fAnalysisManager(aAnalysisManager)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayTelSteppingAction::~XrayTelSteppingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayTelSteppingAction::UserSteppingAction(const G4Step*)
{
#ifdef G4ANALYSIS_USE
  if(fAnalysisManager) fAnalysisManager->Step(fpSteppingManager);
#endif
}







