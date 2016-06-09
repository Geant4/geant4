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
/// @file G4MPIbatch.cc
/// @brief MPI batch session

#include "G4MPIbatch.hh"
#include "G4MPImanager.hh"
#include "G4UImanager.hh"
#include "G4UIcommandStatus.hh"
#include <vector>

// --------------------------------------------------------------------------
static void Tokenize(const G4String& str, std::vector<G4String>& tokens)
{
  const char* delimiter = " ";

  str_size pos0 = str.find_first_not_of(delimiter);
  str_size pos = str.find_first_of(delimiter, pos0);

  while (pos != G4String::npos || pos0 != G4String::npos) {
    if (str[pos0] == '\"') {
      pos = str.find_first_of("\"", pos0+1);
      if(pos != G4String::npos) pos++;
    }
    if (str[pos0] == '\'') {
      pos = str.find_first_of("\'", pos0+1);
      if(pos != G4String::npos) pos++;
    }

    tokens.push_back(str.substr(pos0, pos-pos0));
    pos0 = str.find_first_not_of(delimiter, pos);
    pos = str.find_first_of(delimiter, pos0);
  }
}

// --------------------------------------------------------------------------
G4MPIbatch::G4MPIbatch(const G4String& fname, G4bool qbatch)
  : G4VMPIsession(), isOpened(false), isBatchMode(qbatch)
{
  if(isMaster) {
    batchStream.open(fname, std::ios::in);
    if(batchStream.fail()) {
      G4cerr << "cannot open a macro file(" << fname << ")."
             << G4endl;
    } else {
      isOpened = true;
    }
  }
}

// --------------------------------------------------------------------------
G4MPIbatch::~G4MPIbatch()
{
  if(isOpened) batchStream.close();
}

// --------------------------------------------------------------------------
G4String G4MPIbatch::ReadCommand()
{
  enum { BUFSIZE = 4096 };
  static char linebuf[BUFSIZE];

  G4String cmdtotal = "";
  G4bool qcontinued = false;
  while(batchStream.good()) {
    batchStream.getline(linebuf, BUFSIZE);

    G4String cmdline(linebuf);

    // TAB-> ' ' conversion
    str_size nb = 0;
    while ((nb= cmdline.find('\t',nb)) != G4String::npos) {
      cmdline.replace(nb, 1, " ");
    }

    // strip
    cmdline = cmdline.strip(G4String::both);

    // skip null line if single line
    if(!qcontinued && cmdline.size() == 0) continue;

    // '#' is treated as echoing something
    if(cmdline(0) == '#') return cmdline;

    // tokenize...
    std::vector<G4String> tokens;
    Tokenize(cmdline, tokens);
    qcontinued = false;
    for (G4int i = 0; i < G4int(tokens.size()); i++) {
      // string after '#" is ignored
      if(tokens[i](0) == '#' ) break;
      // '\' or '_' is treated as continued line.
      if(tokens[i] == '\\' || tokens[i] == '_' ) {
        qcontinued = true;
        // check nothing after line continuation character
        if( i != G4int(tokens.size())-1) {
          G4Exception("G4MPIbatch::ReadCommand", "MPI002", JustWarning,
            "unexpected character after line continuation character");
        }
        break; // stop parsing
      }
      cmdtotal += tokens[i];
      cmdtotal += " ";
    }

    if(qcontinued) continue; // read the next line

    if(cmdtotal.size() != 0) break;
    if(batchStream.eof()) break;
  }

  // strip again
  cmdtotal = cmdtotal.strip(G4String::both);

  // bypass some commands
  cmdtotal = BypassCommand(cmdtotal);

  // finally,
  if(batchStream.eof() && cmdtotal.size()==0) {
    return "exit";
  }

  return cmdtotal;
}

// --------------------------------------------------------------------------
G4UIsession* G4MPIbatch::SessionStart()
{
  if( isMaster && !isOpened ){ // macro file is not found
    g4MPI-> BcastCommand("exit");
    return NULL;
  }

  G4String newCommand = "", scommand; // newCommand is always "" in slaves

  while(1) {
    if (isMaster) newCommand = ReadCommand();
    // broadcast a new G4 command
    scommand = g4MPI-> BcastCommand(newCommand);
    if(scommand == "exit") {
      g4MPI-> WaitBeamOn();
      return 0;
    }

    // just echo something
    if( scommand(0) == '#') {
      if(G4UImanager::GetUIpointer()-> GetVerboseLevel() == 2) {
        G4cout << scommand << G4endl;
      }
      continue;
    }

    G4int rc = ExecCommand(scommand);
    if (rc != fCommandSucceeded) {
      G4cerr << G4endl << "***** Batch is interupted!! *****" << G4endl;
      break;
    }
  }

  return NULL;
}

