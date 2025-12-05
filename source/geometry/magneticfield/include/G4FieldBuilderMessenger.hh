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
// G4FieldBuilderMessenger
//
// Class description:
//
// Messenger class that defines commands for G4FieldBuilder.

// Author: Ivana Hrivnacova (IJCLab, Orsay), 2024
// --------------------------------------------------------------------
#ifndef G4FIELDBUILDERMESSENGER_HH
#define G4FIELDBUILDERMESSENGER_HH

#include "G4UImessenger.hh"
#include "globals.hh"

class G4FieldBuilder;

class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;

/**
 * @brief G4FieldBuilderMessenger is messenger class that defines
 * commands for G4FieldBuilder.
 */

class G4FieldBuilderMessenger : public G4UImessenger
{
  public:

    /**
     * Standard Constructor and Destructor.
     */
    G4FieldBuilderMessenger(G4FieldBuilder* fieldBuilder);
    ~G4FieldBuilderMessenger() override;

    /**
     * Default constructor, copy constructor and assignment operator not allowed.
     */
    G4FieldBuilderMessenger() = delete;
    G4FieldBuilderMessenger(const G4FieldBuilderMessenger&) = delete;
    G4FieldBuilderMessenger& operator=(const G4FieldBuilderMessenger&) = delete;

    /**
     * Applies command to the associated object.
     */
    void SetNewValue(G4UIcommand* command, G4String newValues) override;

  private:

    /** Associated class object. */
    G4FieldBuilder* fFieldBuilder = nullptr;

    /** Associated commands directory. */
    G4UIdirectory* fDirectory = nullptr;

    // Commands data members

    /** Command: fieldType. */
    G4UIcmdWithAnInteger* fVerboseLevelCmd = nullptr; 
};

#endif
