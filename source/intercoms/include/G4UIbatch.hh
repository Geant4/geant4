// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIbatch.hh,v 1.2 1999-10-29 06:06:43 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $id$

#ifndef G4UIbatch_h
#define G4UIbatch_h 1

#include "globals.hh"
#include "G4UIsession.hh"
#include <fstream.h>

class G4UImanager;

// class description
//  This is a concrete class of G4UIsession.
//  This class is constructed by G4UImanager when the user 
// applies "/control/execute macro_file" command.
//  If the user's application runs in pure batch mode with
// only one fixed-name macro file, he/she can construct this
// class object with giving the file name and use it as other
// ordinary G4UIsession concrete class, i.e. invoke SessionStart()
// from his/her main().

class G4UIbatch : public G4UIsession 
{
  public: // with description
      G4UIbatch(G4String fileName,G4UIsession* prevSession=NULL);
      //  Constructor. 
      //  "prevSession" must be NULL if this class is constructed
      // from main().
  public:
      ~G4UIbatch();

      G4UIsession * SessionStart();
      void PauseSessionStart(G4String Prompt);
  
  private:
      G4UImanager * UImanager;
      G4UIsession * previousSession;
      ifstream macroFile;
      G4String macroFileName;
      G4bool openFailed;
};



#endif

