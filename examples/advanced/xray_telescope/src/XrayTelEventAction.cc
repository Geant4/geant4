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
// * MODULE:            XrayTelEventAction.cc                           *     
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
// - First implementation of xray_telescope event action
// - Based on Chandra and XMM models 
//
//
// **********************************************************************

#include "G4EventManager.hh"

#include "XrayTelEventAction.hh"
#include "XrayTelEventActionMessenger.hh"

#ifdef G4ANALYSIS_USE
#include "XrayTelAnalysisManager.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayTelEventAction::XrayTelEventAction(XrayTelAnalysisManager* aAnalysisManager)
:fAnalysisManager(aAnalysisManager)
,fDrawFlag("all")
,fEventMessenger(NULL)
{
  fEventMessenger = new XrayTelEventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayTelEventAction::~XrayTelEventAction()
{
  delete fEventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayTelEventAction::BeginOfEventAction(const G4Event*)
{
#ifdef G4ANALYSIS_USE
  if(fAnalysisManager) fAnalysisManager->BeginOfEvent();
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayTelEventAction::EndOfEventAction(const G4Event*)
{
#ifdef G4ANALYSIS_USE
  if(fAnalysisManager) fAnalysisManager->EndOfEvent(fpEventManager->GetConstCurrentEvent(),fDrawFlag);
#endif
}















































