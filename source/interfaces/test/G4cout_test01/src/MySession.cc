// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MySession.cc,v 1.1 1999-01-08 16:32:09 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------------

#include "MySession.hh"

MySession::MySession() : logFileName("Mylog")
{ 
  G4cout <<"G4cout in MySession constructor" << endl;
     const char* theFileName = logFileName;
     logFile.open((char*)theFileName);
}

MySession::~MySession() {// close file for logging;
     logFile.close();
}

MySession * MySession::SessionStart() { return NULL; }

G4int MySession::ReceiveG4cout(G4String coutString)
{
  logFile << coutString << flush;
  return 0;
}

G4int MySession::ReceiveG4cerr(G4String cerrString)
{
  cerr << cerrString << flush;
  logFile << cerrString << flush;
  return 0;
}                                                                       
