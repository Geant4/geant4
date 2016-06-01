// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIWo.hh,v 2.3 1998/10/23 14:42:27 barrand Exp $
// GEANT4 tag $Name: geant4-00 $
//
#ifndef G4UIWo_h
#define G4UIWo_h 

#if defined(G4UI_BUILD_WO_SESSION) || defined(G4UI_USE_WO)

#include "G4UIsession.hh"

class G4WoMessenger;

class G4UIWo : public G4UIsession {
public:
  G4UIWo(int,char**);
  ~G4UIWo();
  G4UIsession* SessionStart();
  void Prompt(G4String);
  void SessionTerminate();
  void PauseSessionStart(G4String);
  G4int ReceiveG4cout(G4String);
  G4int ReceiveG4cerr(G4String);
private:
  void SecondaryLoop(G4String);
  G4WoMessenger* woMessenger;
};

#endif

#endif

