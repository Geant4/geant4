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
//
//
//---------------------------------------------------------------
//  GEANT 4 class header file
//
//  GFlashShowerModelMessenger
//
//  Class description:
//
//  Messenger for the GFlash parameterisation shower model control.

//
// Author: Joanna Weng - 9.11.04
//---------------------------------------------------------------

#ifndef GFlashShowerModelMessenger_h
#define GFlashShowerModelMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class GFlashShowerModel;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;

class GFlashShowerModelMessenger: public G4UImessenger
{
  public:

    GFlashShowerModelMessenger(GFlashShowerModel * myModel);
    ~GFlashShowerModelMessenger();
  
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);
  
  private:

    GFlashShowerModel* myModel;
    G4UIdirectory*   myParaDir;
    G4UIcmdWithAnInteger*  FlagCmd;
    G4UIcmdWithAnInteger*  ContCmd; // Containment Check
    G4UIcmdWithADouble*   StepInX0Cmd;
    G4UIcmdWithADoubleAndUnit*   EmaxCmd;
    G4UIcmdWithADoubleAndUnit*   EminCmd;
    G4UIcmdWithADoubleAndUnit*   EkillCmd;
};

#endif
