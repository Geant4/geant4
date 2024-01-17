//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#include "Par04ParallelMessenger.hh"
#include <CLHEP/Units/SystemOfUnits.h>   // for pi
#include <G4ApplicationState.hh>         // for G4State_PreInit, G4State_Idle
#include <G4ThreeVector.hh>              // for G4ThreeVector
#include <G4Types.hh>                    // for G4bool, G4double, G4int
#include <G4UIcommand.hh>                // for G4UIcommand
#include <G4UImessenger.hh>              // for G4UImessenger
#include <G4UIparameter.hh>              // for G4UIparameter
#include <istream>                       // for basic_istream, basic_istream...
#include <string>                        // for operator>>
#include "G4UIcmdWithAnInteger.hh"       // for G4UIcmdWithAnInteger
#include "G4UIcmdWithoutParameter.hh"    // for G4UIcmdWithoutParameter
#include "G4UIcmdWithABool.hh"           // for G4UIcmdWithABool
#include "G4UIdirectory.hh"              // for G4UIdirectory
#include "Par04ParallelFullWorld.hh"         // for Par04ParallelFullWorld

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04ParallelMessenger::Par04ParallelMessenger(Par04ParallelFullWorld* aParallel)
  : G4UImessenger()
  , fParallel(aParallel)
{
  fExampleDir = new G4UIdirectory("/Par04/");
  fExampleDir->SetGuidance("UI commands specific to this example");

  fParallelDir = new G4UIdirectory("/Par04/parallel/");
  fParallelDir->SetGuidance("Parallel construction UI commands");

  fPrintCmd = new G4UIcmdWithoutParameter("/Par04/parallel/print", this);
  fPrintCmd->SetGuidance("Print current settings.");

  fNbSlicesCmd = new G4UIcmdWithAnInteger("/Par04/parallel/setNbOfSlices", this);
  fNbSlicesCmd->SetGuidance("Set number of slices.");
  fNbSlicesCmd->SetParameterName("NbSlices", false);
  fNbSlicesCmd->SetRange("NbSlices>0");
  fNbSlicesCmd->AvailableForStates(G4State_PreInit);
  fNbSlicesCmd->SetToBeBroadcasted(false);

  fNbRowsCmd = new G4UIcmdWithAnInteger("/Par04/parallel/setNbOfRows", this);
  fNbRowsCmd->SetGuidance("Set number of rows.");
  fNbRowsCmd->SetParameterName("NbRows", false);
  fNbRowsCmd->SetRange("NbRows>0");
  fNbRowsCmd->AvailableForStates(G4State_PreInit);
  fNbRowsCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04ParallelMessenger::~Par04ParallelMessenger()
{
  delete fPrintCmd;
  delete fNbSlicesCmd;
  delete fNbRowsCmd;
  delete fParallelDir;
  delete fExampleDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04ParallelMessenger::SetNewValue(G4UIcommand* aCommand, G4String aNewValue)
{
  if(aCommand == fPrintCmd)
  {
    fParallel->Print();
  }
  else if(aCommand == fNbSlicesCmd)
  {
    fParallel->SetNbOfSlices(fNbSlicesCmd->GetNewIntValue(aNewValue));
  }
  else if(aCommand == fNbRowsCmd)
  {
    fParallel->SetNbOfRows(fNbRowsCmd->GetNewIntValue(aNewValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String Par04ParallelMessenger::GetCurrentValue(G4UIcommand* aCommand)
{
  G4String cv;

  if(aCommand == fNbSlicesCmd)
  {
    cv = fNbSlicesCmd->ConvertToString(fParallel->GetNbOfSlices());
  }
  else if(aCommand == fNbRowsCmd)
  {
    cv = fNbRowsCmd->ConvertToString(fParallel->GetNbOfRows());
  }
  return cv;
}
