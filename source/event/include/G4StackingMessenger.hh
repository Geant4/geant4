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
// $Id: G4StackingMessenger.hh 66892 2013-01-17 10:57:59Z gunter $
//

#ifndef G4StackingMessenger_h
#define G4StackingMessenger_h 1

#include "G4UImessenger.hh"
class G4StackManager;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;

// class description:
//
//  This is a concrete class of G4UImessenger which handles the commands
// for G4StackManager. It has the following commands:
//   /event/stack/
//   /event/stack/status
//   /event/stack/clear
//   /event/stack/verbose

class G4StackingMessenger: public G4UImessenger
{
  public:
    G4StackingMessenger(G4StackManager* fCont);
    ~G4StackingMessenger();
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    G4StackManager * fContainer;
    G4UIdirectory* stackDir;
    G4UIcmdWithoutParameter* statusCmd;
    G4UIcmdWithAnInteger* clearCmd;
    G4UIcmdWithAnInteger* verboseCmd;
};

#endif

