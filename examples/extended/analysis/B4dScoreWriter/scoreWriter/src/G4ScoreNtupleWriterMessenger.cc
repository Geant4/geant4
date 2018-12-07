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
//
// Author: Ivana Hrivnacova, 30/10/2018  (ivana@ipno.in2p3.fr)
// ---------------------------------------------------------------------

#include "G4ScoreNtupleWriterMessenger.hh"
#include "G4ScoreNtupleWriter.hh"
#include "G4UImanager.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4AnalysisUtilities.hh"


G4ScoreNtupleWriterMessenger::G4ScoreNtupleWriterMessenger(
                                G4ScoreNtupleWriter* scoreNtupleWriter)
 : G4UImessenger(),
   fScoreNtupleWriter(scoreNtupleWriter),
   // fDirectory(nullptr),
   fWriterTypeCmd(nullptr),
   fWriterFileNameCmd(nullptr),
   fWriterVerboseCmd(nullptr)
{
// Create command directory and commands only if they do not exist.
// To be improved: G4UImanager::FindDirectory cannot be  used
// as it is declared private.

  // if ( ! G4UImanager::GetUIpointer()->FindDirectory("/score") ) {
  //   fDirectory = new G4UIdirectory("/score/");
  //   fDirectory->SetGuidance("Interactive scoring commands.");
  // }

  fWriterTypeCmd = new G4UIcmdWithAString("/score/writerType",this);
  fWriterTypeCmd->SetGuidance( "Set the ntuple writer output file type.");
  fWriterTypeCmd->SetParameterName("outputType",false);
  G4String candidates;
  for ( auto candidate : 
       { G4AnalysisOutput::kCsv, 
         G4AnalysisOutput::kRoot, 
         G4AnalysisOutput::kXml, 
         G4AnalysisOutput::kNone } ) {
    candidates += G4Analysis::GetOutputName(candidate);
    candidates += " ";  
  }
  fWriterTypeCmd->SetCandidates(candidates);
  fWriterTypeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fWriterFileNameCmd = new G4UIcmdWithAString("/score/writerFileName",this);
  fWriterFileNameCmd->SetGuidance( "Set the ntuple writer output file name.");
  fWriterFileNameCmd->SetParameterName("outputFileName",false);
  fWriterFileNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fWriterVerboseCmd = new G4UIcmdWithAnInteger("/score/writerVerbose",this);
  fWriterVerboseCmd->SetGuidance("Set the  ntuple writer verbose level.");
  fWriterVerboseCmd->SetParameterName("writerVerbose",false);
  fWriterTypeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

G4ScoreNtupleWriterMessenger::~G4ScoreNtupleWriterMessenger()
{
  delete fWriterTypeCmd;
  delete fWriterFileNameCmd;
  delete fWriterVerboseCmd;
  // delete fDirectory;
}

void G4ScoreNtupleWriterMessenger::SetNewValue(G4UIcommand * command,G4String newVal)
{
  if ( command==fWriterTypeCmd ) { 
    fScoreNtupleWriter->SetOutputType(newVal);
  }
  else if ( command==fWriterFileNameCmd ) { 
    fScoreNtupleWriter->SetFileName(newVal);
  }
  else if ( command==fWriterVerboseCmd ) {
    fScoreNtupleWriter->SetVerboseLevel(fWriterVerboseCmd->GetNewIntValue(newVal)); 
  }
}
