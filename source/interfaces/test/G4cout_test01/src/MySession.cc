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
//
// $Id: MySession.cc,v 1.4 2002-11-09 07:17:09 yhajime Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------------

#include "MySession.hh"

MySession::MySession() : logFileName("Mylog")
{ 
  //  G4UImanager * UI = G4UImanager::GetUIpointer();
  G4cout <<"MySession starts G4cout redirection" << G4endl;
     const char* theFileName = logFileName;
     logFile.open((char*)theFileName);
}

MySession::~MySession() {// close file for logging;
     logFile.close();
  G4cout <<"end of MySession" << G4endl;
  UI->SetCoutDestination(NULL);
}

MySession * MySession::SessionStart() { return NULL; }

G4int MySession::ReceiveG4cout(G4String coutString)
{
  //  G4cout << coutString << G4std::flush;
  logFile << coutString << G4std::flush;
  return 0;
}

G4int MySession::ReceiveG4cerr(G4String cerrString)
{
  G4cerr << cerrString << G4std::flush;
  logFile << cerrString << G4std::flush;
  return 0;
}                                                                       
