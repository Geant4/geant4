// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MySession.cc,v 1.2 1999-12-15 14:50:49 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------------

#include "MySession.hh"

MySession::MySession() : logFileName("Mylog")
{ 
  G4cout <<"G4cout in MySession constructor" << G4endl;
     const char* theFileName = logFileName;
     logFile.open((char*)theFileName);
}

MySession::~MySession() {// close file for logging;
     logFile.close();
}

MySession * MySession::SessionStart() { return NULL; }

G4int MySession::ReceiveG4cout(G4String coutString)
{
  logFile << coutString << G4std::flush;
  return 0;
}

G4int MySession::ReceiveG4cerr(G4String cerrString)
{
  G4cerr << cerrString << G4std::flush;
  logFile << cerrString << G4std::flush;
  return 0;
}                                                                       
