// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIsession.cc,v 2.2 1998/10/13 05:36:06 masayasu Exp $
// GEANT4 tag $Name: geant4-00 $
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
