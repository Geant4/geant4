//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
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
      G4UIdirectory*        directory, *subdir1, *subdir2,
                            *subdir3, *subdir4, *subdir5;
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

