// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIsession.cc,v 1.1 1999-01-07 16:09:28 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------------

#include "G4UIsession.hh"

G4UIsession::G4UIsession() {;}

G4UIsession::~G4UIsession() {;}

G4UIsession * G4UIsession::SessionStart() { return NULL; }

void G4UIsession::PauseSessionStart(G4String Prompt) {;}

G4int G4UIsession::ReceiveG4cout(G4String coutString)
{
  cout <<  coutString << flush;
  return 0;
}

G4int G4UIsession::ReceiveG4cerr(G4String cerrString)
{
  cerr <<  cerrString << flush;
  return 0;
}                                                                       
