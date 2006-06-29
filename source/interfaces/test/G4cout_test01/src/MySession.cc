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
// $Id: MySession.cc,v 1.6 2006-06-29 19:10:50 gunter Exp $
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
  //  G4cout << coutString << std::flush;
  logFile << coutString << std::flush;
  return 0;
}

G4int MySession::ReceiveG4cerr(G4String cerrString)
{
  G4cerr << cerrString << std::flush;
  logFile << cerrString << std::flush;
  return 0;
}                                                                       
