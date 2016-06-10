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
// $Id: G4FileMessenger.cc 66310 2012-12-17 11:56:35Z ihrivnac $

// Author: Ivana Hrivnacova, 18/06/2013 (ivana@ipno.in2p3.fr)

#include "G4FileMessenger.hh"
#include "G4VAnalysisManager.hh"
#include "G4AnalysisUtilities.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

using namespace G4Analysis;

//_____________________________________________________________________________
G4FileMessenger::G4FileMessenger(G4VAnalysisManager* manager)
  : G4UImessenger(),
    fManager(manager),
    fSetFileNameCmd(nullptr),
    fSetHistoDirNameCmd(nullptr),
    fSetNtupleDirNameCmd(nullptr)
{  
  fSetFileNameCmd = G4Analysis::make_unique<G4UIcmdWithAString>("/analysis/setFileName",this);
  fSetFileNameCmd->SetGuidance("Set name for the histograms & ntuple file");
  fSetFileNameCmd->SetParameterName("Filename", false);
  fSetFileNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fSetHistoDirNameCmd = G4Analysis::make_unique<G4UIcmdWithAString>("/analysis/setHistoDirName",this);
  fSetHistoDirNameCmd->SetGuidance("Set name for the histograms directory");
  fSetHistoDirNameCmd->SetParameterName("HistoDirName", false);
  fSetHistoDirNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fSetNtupleDirNameCmd = G4Analysis::make_unique<G4UIcmdWithAString>("/analysis/setNtupleDirName",this);
  fSetNtupleDirNameCmd->SetGuidance("Set name for the ntuple directory");
  fSetNtupleDirNameCmd->SetParameterName("NtupleDirName", false);
  fSetNtupleDirNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//_____________________________________________________________________________
G4FileMessenger::~G4FileMessenger()
{}

//
// public functions

//_____________________________________________________________________________
void G4FileMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
  if ( command == fSetFileNameCmd.get() ) {
    G4cout << "Set file name: " << newValues << G4endl;
    fManager->SetFileName(newValues);
  }  
  else if ( command == fSetHistoDirNameCmd.get() ) {
    fManager->SetHistoDirectoryName(newValues);
  }  
  else if ( command == fSetNtupleDirNameCmd.get() ) {
    fManager->SetNtupleDirectoryName(newValues);
  }  
}  
