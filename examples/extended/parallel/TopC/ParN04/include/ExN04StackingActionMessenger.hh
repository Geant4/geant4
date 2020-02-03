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
/// \file ExN04StackingActionMessenger.hh
/// \brief Definition of the ExN04StackingActionMessenger class
//

#ifndef ExN04StackingActionMessenger_h
#define ExN04StackingActionMessenger_h 1

class ExN04StackingAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;

#include "G4UImessenger.hh"
#include "globals.hh"

class ExN04StackingActionMessenger: public G4UImessenger
{
  public:
    ExN04StackingActionMessenger(ExN04StackingAction* msa);
    ~ExN04StackingActionMessenger();
    
  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    ExN04StackingAction * myAction;
    
  private: //commands
    G4UIcmdWithAnInteger * muonCmd;
    G4UIcmdWithAnInteger * isomuonCmd;
    G4UIcmdWithAnInteger * isoCmd;
    G4UIcmdWithADoubleAndUnit * roiCmd;
    
};

#endif


