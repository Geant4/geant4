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

// Author: Ivana Hrivnacova, 29/10/2021  (ivana@ipno.in2p3.fr)

#include "G4ToolsAnalysisMessenger.hh"
#include "G4ToolsAnalysisManager.hh"
#include "G4AnalysisUtilities.hh"
#include "G4AnalysisMessengerHelper.hh"

#include "G4UIcommand.hh"
#include "G4UIparameter.hh"

using namespace G4Analysis;

//_____________________________________________________________________________
G4ToolsAnalysisMessenger::G4ToolsAnalysisMessenger(G4ToolsAnalysisManager* manager)
  : G4UImessenger(),
    fManager(manager)
{
  // Create get commands
  G4AnalysisMessengerHelper helper("h1");
  fGetH1Cmd = helper.CreateGetCommand(this);
  fGetH1VectorCmd = helper.CreateGetVectorCommand(this);

  helper.SetHnType("h2");
  fGetH2Cmd = helper.CreateGetCommand(this);
  fGetH2VectorCmd = helper.CreateGetVectorCommand(this);

  helper.SetHnType("h3");
  fGetH3Cmd = helper.CreateGetCommand(this);
  fGetH3VectorCmd = helper.CreateGetVectorCommand(this);

  helper.SetHnType("p1");
  fGetP1Cmd = helper.CreateGetCommand(this);
  fGetP1VectorCmd = helper.CreateGetVectorCommand(this);

  helper.SetHnType("p2");
  fGetP2Cmd = helper.CreateGetCommand(this);
  fGetP2VectorCmd = helper.CreateGetVectorCommand(this);
}

//_____________________________________________________________________________
G4ToolsAnalysisMessenger::~G4ToolsAnalysisMessenger() = default;

///
// public functions
//

//_____________________________________________________________________________
G4String G4ToolsAnalysisMessenger::GetCurrentValue (G4UIcommand* command)
{
  if ( command == fGetH1Cmd.get() ) return fH1Value;
  if ( command == fGetH2Cmd.get() ) return fH2Value;
  if ( command == fGetH3Cmd.get() ) return fH3Value;
  if ( command == fGetP1Cmd.get() ) return fP1Value;
  if ( command == fGetP2Cmd.get() ) return fP2Value;

  if ( command == fGetH1VectorCmd.get() ) return fH1VectorValue;
  if ( command == fGetH2VectorCmd.get() ) return fH2VectorValue;
  if ( command == fGetH3VectorCmd.get() ) return fH3VectorValue;
  if ( command == fGetP1VectorCmd.get() ) return fP1VectorValue;
  if ( command == fGetP2VectorCmd.get() ) return fP2VectorValue;

  return "";
}

//_____________________________________________________________________________
void G4ToolsAnalysisMessenger::SetNewValue(G4UIcommand* command, G4String value)
{
  auto id = G4UIcommand::ConvertToInt(value);
  if ( command == fGetH1Cmd.get() ) {
    fH1Value = GetHnAddress(id, fManager->fH1Manager);
  }
  else if ( command == fGetH2Cmd.get() ) {
    fH2Value = GetHnAddress(id, fManager->fH2Manager);
  }
  else if ( command == fGetH3Cmd.get() ) {
    fH3Value = GetHnAddress(id, fManager->fH3Manager);
  }
  else if ( command == fGetP1Cmd.get() ) {
    fP1Value = GetHnAddress(id, fManager->fP1Manager);
  }
  else if ( command == fGetP2Cmd.get() ) {
    fP2Value = GetHnAddress(id, fManager->fP2Manager);
  }
  if ( command == fGetH1VectorCmd.get() ) {
    fH1VectorValue = GetHnVectorAddress(fManager->fH1Manager);
  }
  else if ( command == fGetH2VectorCmd.get() ) {
    fH2VectorValue = GetHnVectorAddress(fManager->fH2Manager);
  }
  else if ( command == fGetH3VectorCmd.get() ) {
    fH3VectorValue = GetHnVectorAddress(fManager->fH3Manager);
  }
  else if ( command == fGetP1VectorCmd.get() ) {
    fP1VectorValue = GetHnVectorAddress(fManager->fP1Manager);
  }
  else if ( command == fGetP2VectorCmd.get() ) {
    fP2VectorValue = GetHnVectorAddress(fManager->fP2Manager);
  }
}
