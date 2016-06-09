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
// $Id: SteppingAction.cc,v 1.6 2004/07/23 15:39:39 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"

#include "G4SteppingManager.hh"
#include "G4VProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(RunAction* run, EventAction* event, HistoManager* histo)
:runAction(run), eventAction(event), histoManager(histo)
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

  G4double charge  = aStep->GetTrack()->GetDefinition()->GetPDGCharge();
  G4double steplen = aStep->GetStepLength();
  if (charge != 0.) histoManager->FillHisto(3,steplen);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


