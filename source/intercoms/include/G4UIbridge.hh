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
// G4UIbridge
//
// Class description:
//
// This class is to be used for MT mode.
// Register a particular thread-local G4UImanager with a UI command
// directory name. When a UI command is issued in the master thread
// that starts with this registered directory name, it is immediately
// forwarded to the registered G4UImanager. Such forwarded command
// is not processed in the master thread nor by other worker thread

// Author: A.Dotti, 2013
// --------------------------------------------------------------------
#ifndef G4UIbridge_hh
#define G4UIbridge_hh 1

#include "globals.hh"

class G4UImanager;

class G4UIbridge
{
  public:

    G4UIbridge(G4UImanager* localUI, G4String dir);
    ~G4UIbridge() = default;

    G4int ApplyCommand(const G4String& aCmd);

    inline G4UImanager* LocalUI() const { return localUImanager; }
    inline const G4String& DirName() const { return dirName; }
    inline G4int DirLength() const { return (G4int)dirName.length(); }

  private:

    G4UImanager* localUImanager = nullptr;
    G4String dirName;
};

#endif
