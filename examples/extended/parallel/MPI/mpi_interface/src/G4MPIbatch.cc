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
// $Id: G4MPIbatch.cc,v 1.1 2007/11/16 14:05:41 kmura Exp $
// $Name: geant4-09-02 $
//
// ====================================================================
//   G4MPIsession.cc
//
//                                         2007 Q
// ====================================================================
#include "G4MPIbatch.hh"
#include "G4MPImanager.hh"
#include "G4UIcommandStatus.hh"

// ====================================================================
//
// class description
//
// ====================================================================

////////////////////////////////////////////////////////////
G4MPIbatch::G4MPIbatch(const G4String& fname, G4bool qbatch)
  : G4VMPIsession(),
    isOpened(false), isBatchMode(qbatch)
////////////////////////////////////////////////////////////
{
  if(isMaster) {
    batchStream.open(fname, std::ios::in);
    if(batchStream.fail()) {
      G4cerr << "cannot open a macro file(" << fname << ")."
             << G4endl;
    } else {
      isOpened= true;
    }
  }
}


/////////////////////////////
G4MPIbatch::~G4MPIbatch()
/////////////////////////////
{
  if(isOpened) batchStream.close();
}


//////////////////////////////////
G4String G4MPIbatch::ReadCommand()
//////////////////////////////////
{
  enum { BUFSIZE= 256 };
  char linebuf[BUFSIZE];
  
  while(batchStream.good()) {
    batchStream.getline(linebuf, BUFSIZE);
    
    G4String cmdline(linebuf);
    cmdline= cmdline.strip(G4String::both);
    cmdline= cmdline.strip(G4String::both, '\011'); // remove TAB
    cmdline= TruncateCommand(cmdline);
    
    str_size ic= cmdline.find_first_of('#');
    if(ic != G4String::npos) {
      cmdline= cmdline(0, ic);
    }
    
    if(batchStream.eof()) {
      if(cmdline.size()==0) {
        return "exit";
      } else {
        return cmdline;
      }
    }

    if(cmdline.size()==0) continue; // skip null line
    return cmdline;    
  }  

  return "exit"; // dummy
}


///////////////////////////////////////
G4UIsession* G4MPIbatch::SessionStart()
///////////////////////////////////////
{
  G4String newCommand="", scommand; // newCommand is always "" in slaves

  if(isMaster) newCommand= ReadCommand();
  // broadcast a new G4 command
  scommand= g4MPI-> BcastCommand(newCommand); 
  if(scommand == "exit" ) return 0;
  
  while(1){
    G4int rc= ExecCommand(scommand);
    if(rc != fCommandSucceeded) break;
    
    if(isMaster) newCommand= ReadCommand();
    scommand= g4MPI-> BcastCommand(newCommand);
    if(scommand == "exit" ) {
      if(isBatchMode) {
        if(g4MPI-> CheckThreadStatus()) {
          g4MPI-> JoinBeamOnThread();
          break;
        }
      } else {
        break;
      }
    }
  }
  
  return 0;
}

