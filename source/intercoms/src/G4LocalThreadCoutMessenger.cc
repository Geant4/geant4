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
// G4LocalThreadCoutMessenger
//
// Author: M.Asai, 2013
// --------------------------------------------------------------------

#include "G4LocalThreadCoutMessenger.hh"

#include "G4UImanager.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIparameter.hh"
#include "G4Tokenizer.hh"

// --------------------------------------------------------------------
G4LocalThreadCoutMessenger::G4LocalThreadCoutMessenger()
{
  coutDir = new G4UIdirectory("/control/cout/");
  coutDir->SetGuidance("Control cout/cerr for local thread.");

  coutFileNameCmd = new G4UIcommand("/control/cout/setCoutFile", this);
  coutFileNameCmd->SetGuidance(
    "Send G4cout stream to a file dedicated to a thread. ");
  coutFileNameCmd->SetGuidance(
    "To have a display output, use special keyword \"**Screen**\".");
  coutFileNameCmd->SetGuidance(
    "If append flag is true output is appended to file,");
  coutFileNameCmd->SetGuidance("otherwise file output is overwritten.");
  coutFileNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  auto* pp = new G4UIparameter("fileName", 's', true);
  pp->SetDefaultValue("**Screen**");
  coutFileNameCmd->SetParameter(pp);
  pp = new G4UIparameter("append", 'b', true);
  pp->SetDefaultValue(1);
  coutFileNameCmd->SetParameter(pp);

  cerrFileNameCmd = new G4UIcommand("/control/cout/setCerrFile", this);
  cerrFileNameCmd->SetGuidance(
    "Send G4cerr stream to a file dedicated to a thread. ");
  cerrFileNameCmd->SetGuidance(
    "To have a display output, use special keyword \"**Screen**\".");
  cerrFileNameCmd->SetGuidance(
    "If append flag is true output is appended to file,");
  cerrFileNameCmd->SetGuidance("otherwise file output is overwritten.");
  cerrFileNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  pp = new G4UIparameter("fileName", 's', true);
  pp->SetDefaultValue("**Screen**");
  cerrFileNameCmd->SetParameter(pp);
  pp = new G4UIparameter("append", 'b', true);
  pp->SetDefaultValue(1);
  cerrFileNameCmd->SetParameter(pp);

  bufferCoutCmd = new G4UIcmdWithABool("/control/cout/useBuffer", this);
  bufferCoutCmd->SetGuidance("Send cout and/or cerr stream to a buffer.");
  bufferCoutCmd->SetGuidance(
    "The buffered text will be printed at the end of the job");
  bufferCoutCmd->SetGuidance(
    "for each thread at a time, so that output of each thread is grouped.");
  bufferCoutCmd->SetGuidance(
    "This command has no effect if output goes to a file.");
  bufferCoutCmd->SetParameterName("flag", true);
  bufferCoutCmd->SetDefaultValue(true);
  bufferCoutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  prefixCmd = new G4UIcmdWithAString("/control/cout/prefixString", this);
  prefixCmd->SetGuidance(
    "Set the prefix string for each cout/cerr line from a thread.");
  prefixCmd->SetParameterName("prefix", true);
  prefixCmd->SetDefaultValue("G4WT");
  prefixCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  ignoreCmd =
    new G4UIcmdWithAnInteger("/control/cout/ignoreThreadsExcept", this);
  ignoreCmd->SetGuidance("Omit cout from threads except the specified one.");
  ignoreCmd->SetGuidance("This command takes effect only if cout destination "
                         "is screen without buffering.");
  ignoreCmd->SetGuidance(
    "If specified thread ID is greater than the number of threads,");
  ignoreCmd->SetGuidance(
    "no cout is displayed from worker threads. -1 to reset.");
  ignoreCmd->SetGuidance("This command does not affect to cerr.");
  ignoreCmd->SetParameterName("threadID", true);
  ignoreCmd->SetDefaultValue(0);
  ignoreCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  ignoreInitCmd =
    new G4UIcmdWithABool("/control/cout/ignoreInitializationCout", this);
  ignoreInitCmd->SetGuidance("Omit cout from threads during initialization, as "
                             "they should be identical to the master thread.");
  ignoreInitCmd->SetGuidance("This command takes effect only if cout "
                             "destination is screen without buffering.");
  ignoreInitCmd->SetGuidance("This command does not affect to cerr.");
  ignoreInitCmd->SetParameterName("IgnoreInit", true);
  ignoreInitCmd->SetDefaultValue(true);
  ignoreInitCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

// --------------------------------------------------------------------
G4LocalThreadCoutMessenger::~G4LocalThreadCoutMessenger()
{
  delete coutFileNameCmd;
  delete cerrFileNameCmd;
  delete bufferCoutCmd;
  delete prefixCmd;
  delete ignoreCmd;
  delete ignoreInitCmd;
  delete coutDir;
}

// --------------------------------------------------------------------
void G4LocalThreadCoutMessenger::SetNewValue(G4UIcommand* command,
                                             G4String newVal)
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(command == coutFileNameCmd)
  {
    G4Tokenizer next(newVal);
    G4String fn = next();
    G4bool af   = StoB(next());
    UI->SetCoutFileName(fn, af);
  }
  else if(command == cerrFileNameCmd)
  {
    G4Tokenizer next(newVal);
    G4String fn = next();
    G4bool af   = StoB(next());
    UI->SetCerrFileName(fn, af);
  }
  else if(command == bufferCoutCmd)
  {
    UI->SetThreadUseBuffer(StoB(newVal));
  }
  else if(command == prefixCmd)
  {
    UI->SetThreadPrefixString(newVal);
  }
  else if(command == ignoreCmd)
  {
    UI->SetThreadIgnore(StoI(newVal));
  }
  else if(command == ignoreInitCmd)
  {
    UI->SetThreadIgnoreInit(StoB(newVal));
  }
}
