// $Id: G4PersistencyCenterMessenger.hh,v 1.1 2002-11-24 13:45:23 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4PersistencyCenterMessenger.hh
//
// History:
//   01.07.18  Youhei Morita  Initial creation (with "fadsclass")

#ifndef PERSISTENCY_CENTER_MESSENGER_HH
#define PERSISTENCY_CENTER_MESSENGER_HH 1

#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4PersistencyCenter.hh"

// Class inherited:
#include "G4UImessenger.hh"

// Class Description:
//   User interface messenger class to interface G4PersistencyCenter

class G4PersistencyCenterMessenger
 : public G4UImessenger
{
    public: // With description
      G4PersistencyCenterMessenger(G4PersistencyCenter* p);
      // Constructor

      ~G4PersistencyCenterMessenger();
      // Destructor

    public: // With description
      void SetNewValue(G4UIcommand* command, G4String newValues);
      // User interface for setting a new value

      G4String GetCurrentValue(G4UIcommand* command);
      // User interface for getting a value

    private:
      std::string PopWord(std::string text, int n, std::string delim);
      // Parse text and returns the n-th words separated by delim

    private:
      G4PersistencyCenter*    pc;
      G4UIdirectory*        directory;
      G4UIcmdWithAnInteger* verboseCmd;
      G4UIcmdWithAString*   select;
      G4UIcmdWithAString*   regHitIO;
      std::vector<std::string>         wrObj;
      std::vector<std::string>         rdObj;
      std::vector<G4UIcmdWithAString*> storeObj;
      std::vector<G4UIcmdWithAString*> setWrFile;
      std::vector<G4UIcmdWithAString*> setRdFile;
      G4UIcmdWithoutParameter*         printAll;

}; // End of class G4PersistencyCenterMessenger

#endif

