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
//

#ifndef G4InteractorMessenger_h
#define G4InteractorMessenger_h 1

#include "G4UImessenger.hh"

class G4VInteractiveSession;
class G4UIdirectory;
class G4UIcommand;

class G4InteractorMessenger : public G4UImessenger 
{
public:
  G4InteractorMessenger (G4VInteractiveSession* session);
  virtual ~G4InteractorMessenger ();
  void SetNewValue(G4UIcommand* command,G4String newValue);
private:
  G4VInteractiveSession* session;
  G4UIdirectory* interactorDirectory;
  G4UIcommand* addMenu;
  G4UIcommand* addButton;
  G4UIcommand* addIcon;
  G4UIcommand* defaultIcons;
  G4UIcommand* sys;
  G4UIcommand* outputStyle;
};

#endif
