// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIsession.hh,v 1.3 1999-12-15 14:50:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $id$

#ifndef G4UIsession_h
#define G4UIsession_h 1

#include "G4coutDestination.hh"
#include "globals.hh"

// class description:
//
//  This is a base class of all (G)UI session.
//  SessionStart() method should be called to start the session.
//

class G4UIsession : public G4coutDestination
{
  // Base class of UI/GUI session
  
  public:
      G4UIsession();
      virtual ~G4UIsession();

      virtual G4UIsession * SessionStart();
      // This method will be invoked by main().
      // Optionally, it can be invoked by another session.
      
      virtual void PauseSessionStart(G4String Prompt);
      // This method will be invoked by G4UImanager
      // when G4kernel becomes to Pause state.
      
      virtual G4int ReceiveG4cout(G4String coutString);
      virtual G4int ReceiveG4cerr(G4String cerrString);
      // These two methods will be invoked by G4strstreambuf.

};



#endif

