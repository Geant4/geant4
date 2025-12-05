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
// G4FieldSetupMessenger
//
// Class description:
//
// Messenger class that defines commands for G4FieldSetup.
//
// Implements commands:
// - /field/update

// Author: Ivana Hrivnacova (IJCLab, Orsay), 2024
// --------------------------------------------------------------------
#ifndef G4FIELDSETUPMESSENGER_HH
#define G4FIELDSETUPMESSENGER_HH

#include "G4UImessenger.hh"
#include "globals.hh"

class G4FieldSetup;

class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;

/**
 * @brief G4FieldSetupMessenger is a messenger class that defines
 * commands for G4FieldSetup.
 */

class G4FieldSetupMessenger : public G4UImessenger
{
  public:

    /**
     * Standard Constructor and Destructor.
     */
    G4FieldSetupMessenger(G4FieldSetup* fieldSetup);
    ~G4FieldSetupMessenger() override;

    /**
     * Default constructor, copy constructor and assignment operator not allowed.
     */
    G4FieldSetupMessenger() = delete;
    G4FieldSetupMessenger(const G4FieldSetupMessenger&) = delete;
    G4FieldSetupMessenger& operator=(const G4FieldSetupMessenger&) = delete;

    /**
     * Applies command to the associated object.
     */
    void SetNewValue(G4UIcommand* command, G4String newValues) override;

  private:

    /** Associated class object. */
    G4FieldSetup* fFieldSetup = nullptr;

    // Commands data members

    /** Command: update. */
    G4UIcmdWithoutParameter* fUpdateCmd = nullptr; 
};

#endif
