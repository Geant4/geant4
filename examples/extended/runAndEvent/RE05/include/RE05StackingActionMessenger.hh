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
/// \file RE05/include/RE05StackingActionMessenger.hh
/// \brief Definition of the RE05StackingActionMessenger class
//

#ifndef RE05StackingActionMessenger_h
#define RE05StackingActionMessenger_h 1

class RE05StackingAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;

#include "G4UImessenger.hh"
#include "globals.hh"

class RE05StackingActionMessenger: public G4UImessenger
{
  public:
    RE05StackingActionMessenger(RE05StackingAction* msa);
    virtual ~RE05StackingActionMessenger();
    
  public:
    virtual void SetNewValue(G4UIcommand * command,G4String newValues);
    virtual G4String GetCurrentValue(G4UIcommand * command);

  private:
    RE05StackingAction * fMyAction;
    
    //commands
    G4UIcmdWithAnInteger* fMuonCmd;
    G4UIcmdWithAnInteger* fIsoMuonCmd;
    G4UIcmdWithAnInteger* fIsoCmd;
    G4UIcmdWithADoubleAndUnit* fRoiCmd;  
};

#endif
