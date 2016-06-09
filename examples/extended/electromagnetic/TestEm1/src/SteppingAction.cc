//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: SteppingAction.cc,v 1.3 2003/11/07 15:38:29 maire Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"

#include "G4SteppingManager.hh"
#include "G4VProcess.hh"

#ifdef G4ANALYSIS_USE
 #include "AIDA/IHistogram1D.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(RunAction* RuAct, EventAction* EvAct)
  :runAction(RuAct),eventAction(EvAct)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4double EdepStep = aStep->GetTotalEnergyDeposit();
  if (EdepStep > 0.) {  runAction->AddEdep(EdepStep);
                      eventAction->AddEdep(EdepStep);
  }
  const G4VProcess* process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
  if (process) runAction->CountProcesses(process->GetProcessName());


#ifdef G4ANALYSIS_USE
  G4double charge  = aStep->GetTrack()->GetDefinition()->GetPDGCharge();
  G4double steplen = aStep->GetStepLength();
  if (charge != 0.) runAction->GetHisto(2)->fill(steplen);
#endif                    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


