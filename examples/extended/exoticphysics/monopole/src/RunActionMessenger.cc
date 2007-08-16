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
// $Id: RunActionMessenger.cc,v 1.1 2007-08-16 10:32:04 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunActionMessenger.hh"

#include "RunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunActionMessenger::RunActionMessenger(RunAction * ra): runaction(ra)
{   
  actDir = new G4UIdirectory("/testex/run/");
  actDir->SetGuidance("run commands");
      
  binSizeCmd = new G4UIcmdWithADouble("/testex/run/binSize", this);
  binSizeCmd->SetGuidance("Set bin size.");
  binSizeCmd->SetParameterName("binSize", false);
  binSizeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  verboseCmd = new G4UIcmdWithAnInteger("/testex/run/verbose", this);
  verboseCmd->SetGuidance("Set verbose level.");
  verboseCmd->SetParameterName("verbose", false);
  verboseCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	histoNameCmd = new G4UIcmdWithAString("/testex/run/HistoName", this);
  histoNameCmd->SetGuidance("Set histo file name.");
  histoNameCmd->SetParameterName("histoName", false);
  histoNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	histoTypeCmd = new G4UIcmdWithAString("/testex/run/HistoType", this);
  histoTypeCmd->SetGuidance("Set histo file type.");
  histoTypeCmd->SetParameterName("histoType", false);
  histoTypeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RunActionMessenger::~RunActionMessenger()
{
  delete binSizeCmd;
  delete verboseCmd;
	delete histoNameCmd;
	delete histoTypeCmd;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  if( command == binSizeCmd )
   { runaction->SetBinSize(binSizeCmd->GetNewDoubleValue(newValue));}

  if( command == verboseCmd )
   { runaction->SetVerbose(verboseCmd->GetNewIntValue(newValue));}

  if( command == histoNameCmd )
   { runaction->SetHistoName(newValue);}

  if( command == histoTypeCmd )
   { runaction->SetHistoType(newValue);}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
