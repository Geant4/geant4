// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MySession.hh,v 1.2 1999-12-15 14:50:48 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $id$

#ifndef MySession_h
#define MySession_h 1

#include "globals.hh"
#include "G4UIsession.hh"
#include "g4std/fstream"

class MySession : public G4UIsession 
{
  // Base class of UI/GUI session
  
  public:
      MySession();
      ~MySession();

      MySession * SessionStart();
      // This method will be invoked by main().
      // Optionally, it can be invoked by another session.
      
      //void PauseSessionStart(G4String Prompt);
      // This method will be invoked by G4UImanager
      // when G4kernel becomes to Pause state.
      
      G4int ReceiveG4cout(G4String coutString);
      G4int ReceiveG4cerr(G4String cerrString);
      // These two methods will be invoked by G4strstreambuf.
      
      G4std::ofstream logFile;
      G4String logFileName;

};



#endif

