// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIbatch.hh,v 1.1 1999-01-07 16:09:20 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $id$

#ifndef G4UIbatch_h
#define G4UIbatch_h 1

#include "globals.hh"
#include "G4UIsession.hh"
#include <fstream.h>

class G4UImanager;

class G4UIbatch : public G4UIsession 
{
  public:
      G4UIbatch(G4String fileName,G4UIsession* prevSession=NULL);
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

