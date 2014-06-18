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
/// @file G4VMPIsession.cc
/// @brief A base class for MPI sessions

#include "mpi.h"
#include "G4UIcommand.hh"
#include "G4UImanager.hh"
#include "G4MPImanager.hh"
#include "G4VMPIsession.hh"

// --------------------------------------------------------------------------
G4VMPIsession::G4VMPIsession()
  : G4VBasicShell()
{
  g4mpi_ = G4MPImanager::GetManager();

  is_master_ = g4mpi_-> IsMaster();
  is_slave_  = g4mpi_-> IsSlave();
  rank_ = g4mpi_-> GetRank();

}

// --------------------------------------------------------------------------
G4VMPIsession::~G4VMPIsession()
{
}

// --------------------------------------------------------------------------
G4bool G4VMPIsession::GetHelpChoice(G4int& aval)
{
  G4cin >> aval;
  if( !G4cin.good() ){
    G4cin.clear();
    G4cin.ignore(30,'\n');
    return false;
  }
  return true;
}

// --------------------------------------------------------------------------
void G4VMPIsession::ExitHelp() const
{
  char temp[100];
  G4cin.getline(temp, 100);
}

// --------------------------------------------------------------------------
void G4VMPIsession::PauseSessionStart(const G4String&)
{
}

// --------------------------------------------------------------------------
G4int G4VMPIsession::ExecCommand(const G4String& acommand)
{
  if( acommand.length() < 2 ) return fCommandSucceeded;

  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4int returnVal = 0;

  G4String command = BypassCommand(acommand);

  // "/mpi/beamOn is threaded out.
  if( command(0,11) == "/mpi/beamOn" ) {
    g4mpi_-> ExecuteBeamOnThread(command);
    returnVal = fCommandSucceeded;
  } else if( command(0,12) == "/mpi/.beamOn" ) { // care for beamOn
    G4bool threadStatus = g4mpi_-> CheckThreadStatus();
    if ( threadStatus ) { // still /run/beamOn is active
      if( is_master_ ) {
        G4cout << "G4MPIsession:: beamOn is still running." << G4endl;
      }
      returnVal = fCommandSucceeded;
    } else {
      returnVal = UI-> ApplyCommand(command);
    }
  } else { // normal command
    returnVal = UI-> ApplyCommand(command);
  }

  G4int paramIndex = returnVal % 100;
  // 0 - 98 : paramIndex-th parameter is invalid
  // 99     : convination of parameters is invalid
  G4int commandStatus = returnVal - paramIndex;

  G4UIcommand* cmd = 0;
  if( commandStatus != fCommandSucceeded ) {
    cmd = FindCommand(command);
  }

  switch( commandStatus ) {
  case fCommandSucceeded:
    break;
  case fCommandNotFound:
    G4cerr << "command <" << UI-> SolveAlias(command)
           << "> not found" << G4endl;
    break;
  case fIllegalApplicationState:
    G4cerr << "illegal application state -- command refused" << G4endl;
    break;
  case fParameterOutOfRange:
    // ...
    break;
  case fParameterOutOfCandidates:
    G4cerr << "Parameter is out of candidate list (index "
           << paramIndex << ")" << G4endl;
    G4cerr << "Candidates : "
           << cmd->GetParameter(paramIndex)-> GetParameterCandidates()
           << G4endl;
    break;
  case fParameterUnreadable:
    G4cerr << "Parameter is wrong type and/or is not omittable (index "
           << paramIndex << ")" << G4endl;
    break;
  case fAliasNotFound:
    // ...
    break;
  default:
    G4cerr << "command refused (" << commandStatus << ")" << G4endl;
  }

  return returnVal;
}

// --------------------------------------------------------------------------
G4String G4VMPIsession::TruncateCommand(const G4String& command) const
{
  // replace "//" with "/" in G4command
  G4String acommand = command;
  G4String strarg;

  str_size iarg = acommand.find(' ');
  if( iarg != G4String::npos ) {
    strarg = acommand(iarg, acommand.size()-iarg);
    acommand = acommand(0,iarg);
  }

  str_size idx;
  while( (idx = acommand.find("//")) != G4String::npos)  {
    G4String command1 = acommand(0,idx+1);
    G4String command2 = acommand(idx+2, acommand.size()-idx-2);
    acommand = command1 + command2;
  }

  acommand += strarg;

  return acommand;
}

// --------------------------------------------------------------------------
G4String G4VMPIsession::BypassCommand(const G4String& command) const
{
  // bypass some commands
  // * /mpi/beamOn
  //    -> /mpi/.beamOn (batch session)
  //    -> /mpi/.beamOn (MT mode)
  //
  // * /run/beamOn
  //    -> /mpi/.beamOn (batch session)
  //    -> /mpi/beamOn  (interactive session)
  //    -> /mpi/.beamOn (MT mode)
  //
  // * /control/execute -> /mpi/execute

  G4String acommand = command;

  // /mpi/beamOn
  if( acommand(0,11) == "/mpi/beamOn" ) {
#ifdef G4MULTITHREADED
    acommand = "/mpi/.beamOn";
    if(command.length() > 11) {
      acommand += command.substr(11);
    }
#else
    if( g4mpi_-> IsBatchMode()) {
      acommand = "/mpi/.beamOn";
      if(command.length() > 11) {
        acommand += command.substr(11);
      }
    }
#endif
  }

  // /run/beamOn
  if( acommand(0,11) == "/run/beamOn" ) {
    G4String strarg = "";
    G4bool qget = false;
    G4bool qdone = false;

    for ( str_size idx = 10; idx < command.size(); idx++ ) {
      if( command[idx] == ' ' || command[idx] == '\011' ) {
        qget = true;
        if(qdone) break;
        continue;
      }
      if( qget ) {
        strarg += command[idx];
        qdone = true;
      }
    }

    if( g4mpi_-> IsBatchMode() ) { // batch session
      acommand = "/mpi/.beamOn ";
      if( command.length() > 11 ) acommand += strarg;
    } else { // interactive session
#ifdef G4MULTITHREADED
      if( g4mpi_-> GetVerbose()>0 && is_master_ ) {
        G4cout << "/run/beamOn is overridden by /mpi/.beamOn" << G4endl;
      }
      acommand = "/mpi/.beamOn ";
      if( command.length() > 11 ) acommand += strarg;
#else
      if( g4mpi_-> GetVerbose()>0 && is_master_ ) {
        G4cout << "/run/beamOn is overridden by /mpi/beamOn" << G4endl;
      }
      acommand = "/mpi/beamOn ";
      if( command.length() > 11 ) acommand += strarg;
#endif
    }
  }

  // /control/execute
  if( acommand(0,16) == "/control/execute" ) {
    if( g4mpi_-> GetVerbose()>0 && is_master_ ) {
      G4cout << "/control/execute is overridden by /mpi/execute"
             << G4endl;
    }
    acommand.replace(0, 16, "/mpi/execute    ");
  }

  return acommand;
}

// --------------------------------------------------------------------------
G4int G4VMPIsession::ReceiveG4cout(const G4String& coutString)
{
  g4mpi_-> Print(coutString);
  return 0;
}

// --------------------------------------------------------------------------
G4int G4VMPIsession::ReceiveG4cerr(const G4String& cerrString)
{
  g4mpi_-> Print(cerrString);
  return 0;
}
