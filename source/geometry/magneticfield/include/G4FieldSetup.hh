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
// G4FieldSetup
//
// Class description:
//
// The class for constructing magnetic, electromagnetic and gravity
// fields which strength is defined via G4Field.
//
// The equation of motion of a particle in a field and the
// integration method are set according to the selection in
// G4FieldParameters, as well as other accuracy parameters.
// The default values in G4FieldParameters correspond to defaults
// set in Geant4.

// Author: Ivana Hrivnacova (IJClab, Orsay), 2024.
// --------------------------------------------------------------------
#ifndef G4FIELDSETUP_HH
#define G4FIELDSETUP_HH

#include "G4FieldParameters.hh"
#include "globals.hh"

class G4Field;
class G4FieldParameters;
class G4FieldSetupMessenger;

class G4ChordFinder;
class G4EquationOfMotion;
class G4FieldManager;
class G4MagIntegratorStepper;
class G4LogicalVolume;
class G4VIntegrationDriver;

/**
 * @brief G4FieldSetup is a class for constructing magnetic, electromagnetic
 * and gravity fields which strength is defined via G4Field.
 * The equation of motion of a particle in a field and the integration method
 * are set according to the selection in G4FieldParameters, as well as other
 * accuracy parameters.
 */

class G4FieldSetup
{
  public:

    /**
     * Standard constructor for G4FieldSetup.
     *  @param[in] parameters The field parameters.
     *  @param[in] field Pointer to the field object.
     *  @param[in] lv Optional logical volume where field applies; if
     *             null, global field applies.
     */
    G4FieldSetup(const G4FieldParameters& parameters,
                       G4Field* field,
                       G4LogicalVolume* lv = nullptr);

    /**
     * Default Destructor.
     */
    ~G4FieldSetup();

    /**
     * Default constructor, copy constructor and assignment operator not allowed.
     */
    G4FieldSetup() = delete;
    G4FieldSetup(const G4FieldSetup& right) = delete;
    G4FieldSetup& operator=(const G4FieldSetup& right) = delete;

    /**
     * Clears previously created setup.
     */
    void Clear();

    /**
     * Updates the field setup with new field parameters.
     */
    void Update();

    /**
     * Prints information.
     *  @param[in] verboseLevel Verbosity level; if greater than 1, parameters
     *             are also printed out to standard output.
     *  @param[in] about Optional string.
     */
    void PrintInfo(G4int verboseLevel, const G4String& about = "created");

    /**
     * Setter for the field object.
     */
    inline void SetG4Field(G4Field* field) { fG4Field = field; }

    /**
     * Accessors.
     */
    inline G4Field* GetG4Field() const { return fG4Field; }
    inline G4LogicalVolume* GetLogicalVolume() const { return fLogicalVolume; }
    inline G4EquationOfMotion* GetEquation() const { return fEquation; }
    inline G4MagIntegratorStepper* GetStepper() const { return fStepper; }

  private:

    /**
     * Creates cached magnetic field if const distance is set greater than zero.
     *  @param[in] parameters The field parameters.
     *  @param[in] field Pointer to the field in input.
     *  @returns The pointer to the cached field or the input field otherwise.
     */
    G4Field* CreateCachedField( const G4FieldParameters& parameters,
                                      G4Field* field);

    /**
     * Creates and sets the equation of motion of a particle in a field.
     *  @param[in] equation The equation type.
     *  @returns The pointer to the created equation of motion.
     */
    G4EquationOfMotion* CreateEquation(G4EquationType equation);

    /**
     * Creates and sets the field integration stepper.
     *  @param[in] equation Pointer to the equation of motion.
     *  @param[in] stepper The stepper type.
     *  @returns The pointer to the created integration stepper.
     */
    G4MagIntegratorStepper* CreateStepper(G4EquationOfMotion* equation,
                                          G4StepperType stepper);

    /**
     * Creates and sets the FSAL field integration driver.
     *  @param[in] equation Pointer to the equation of motion.
     *  @param[in] stepper The stepper type.
     *  @param[in] minStep The minimum allowed step.
     *  @returns The pointer to the created FSAL integration driver.
     */
    G4VIntegrationDriver*
    CreateFSALStepperAndDriver(G4EquationOfMotion* equation,
                               G4StepperType stepper, G4double minStep);

    // Methods to update field setup step by step

    /**
     * Creates cached field (if ConstDistance is set).
     */
    void CreateCachedField();

    /**
     * Creates the stepper.
     */
    void CreateStepper();

    /**
     * Creates the chord finder.
     */
    void CreateChordFinder();

    /**
     * Updates the field manager.
     */
    void UpdateFieldManager();

  private:  // data members

    /** Messenger for this class. */
    G4FieldSetupMessenger* fMessenger = nullptr;

    /** Field parameters. */
    const G4FieldParameters& fParameters;

    /** The field manager. */
    G4FieldManager* fFieldManager = nullptr;

    /** The field class object. */
    G4Field* fG4Field = nullptr;

    /** The associated volume (if local field). */
    G4LogicalVolume* fLogicalVolume = nullptr;

    /** The equation of motion. */
    G4EquationOfMotion* fEquation = nullptr;

    /** The magnetic integrator stepper. */
    G4MagIntegratorStepper* fStepper = nullptr;

    /** The magnetic integrator driver. */
    G4VIntegrationDriver* fDriver = nullptr;

    /** Chord finder. */
    G4ChordFinder* fChordFinder = nullptr;
};

#endif
