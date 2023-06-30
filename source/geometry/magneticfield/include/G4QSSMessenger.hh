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
// G4QSSMessenger
//
// Messenger for QSS Integrator driver

// Author: Leandro Gomez Vidal - October 2021
// --------------------------------------------------------------------
#ifndef GEANT4_G4QSSMessenger_H
#define GEANT4_G4QSSMessenger_H

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIdirectory.hh"
#include "G4UImessenger.hh"

class G4QSSMessenger : public G4UImessenger
{
   public:

     G4QSSMessenger();
    ~G4QSSMessenger() override;

     void SetNewValue(G4UIcommand* command, G4String newValues) override;

     /* Hacky - much easier to access G4QSSMessenger from G4QSSDriver than the other way around
      * Multithreading seems to cause weird stuff with Driver/Stepper instances, maybe making
      * some thread-local copies or something, outside of construction? */

     static G4QSSMessenger* instance();

     enum StepperSelection
     {
       None = 0,
       TemplatedDoPri,
       OldRK45,
       G4QSS2
     };
     void selectStepper(const std::string&);
     StepperSelection selectedStepper();

  public:

    G4double dQMin = 0;
    G4double dQRel = 0;
    G4double trialProposedStepModifier = 1.0;

  private:

    StepperSelection _selectedStepper;
    G4UIdirectory* qssCmdDir;
    G4UIcmdWithADoubleAndUnit* dQMinCmd;
    G4UIcmdWithADouble* dQRelCmd;
    G4UIcmdWithAString* stepperSelectorCmd;
    G4UIcmdWithADouble* trialProposedStepModifierCmd;
};

#endif  // GEANT4_G4QSSMessenger_H
