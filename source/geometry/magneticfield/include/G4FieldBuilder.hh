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
// G4FieldBuilder
//
// Class description:
//
// The manager class for building magnetic or other fields
// using the configuration in field parameters.
//
// Purpose: Provide a single 'place' to configure field & integration
//
// - It can configure a global field, and field(s) local to a (logical) volume
// - The parameter values can be configured by the user (else use a default)
// - They can be set/changed via a messenger provided or in the code
//      of the user detector construciton
// - It retains ownership of the following object(s):
//      field parameters and field setups, field
//
// Note MT: an object of the builder class should be created on master only
//      (in DetectorConstruction constructor)
//      The functions SetGlobal/LocalField and ConstructFieldSetup should
//      be called on workers (in DetectorConstruction::ConstructSDandField )
//
// This design/implementation covers the most common use cases.
// It cannot be used to create some complex setups such as
//    - equations templated on the field type,
//    - steppers/drivers templated on the equation and field types.

// Author: Ivana Hrivnacova (IJCLab, Orsay), 2024.
// --------------------------------------------------------------------
#ifndef G4FIELDBUILDER_HH
#define G4FIELDBUILDER_HH

#include "G4Cache.hh"
#include "G4FieldParameters.hh"
#include "globals.hh"

#include <vector>

class G4Field;
class G4FieldBuilderMessenger;
class G4FieldSetup;
class G4LogicalVolume;
class G4EquationOfMotion;
class G4MagIntegratorStepper;

/**
 * @brief G4FieldBuilder is a singleton manager class for building magnetic
 * or other fields, using the configuration in G4FieldParameters.
 */

class G4FieldBuilder
{
  public:

    /**
     * Destructor. Deletes associated messenger.
     */
    ~G4FieldBuilder();

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4FieldBuilder(const G4FieldBuilder& right) = delete;
    G4FieldBuilder& operator=(const G4FieldBuilder& right) = delete;

    // Static access methods

    /**
     * Creates the class instance, if it does not exist; simply returns it
     * on the next calls.
     *  @returns A pointer to the singleton instance.
     */
    static G4FieldBuilder* Instance();

    /**
     * Tells if the singleton instance exists.
     *  @returns true if the singleton instance exists.
     */
    static G4bool IsInstance();

    // Functions for constructing field setup

    /**
     * Creates the local magnetic field parameters (configuration) which can
     * be then configured by the user via UI commands.
     * The parameters are used in geometry only if a local magnetic field is
     * associated with the volumes with the given name.
     *  @param[in] fieldVolName Volume name.
     *  @returns A pointer to the field parameters.
     */
    G4FieldParameters* CreateFieldParameters(const G4String& fieldVolName);

    /**
     * Constructs the setup for all registered fields.
     */
    void ConstructFieldSetup();

    /**
     * Updates the magnetic field. It must be called if the field parameters
     * were changed in other than PreInit> phase.
     */
    void UpdateField();

    /**
     * Reinitialises if geometry has been modified. This method is called
     * by G4RunManager during ReinitializeGeometry().
     */
    void Reinitialize();

    // Set methods

    /**
     * The default field type is set to "kMagnetic". This method should be
     * called for other than magnetic field, in order to update the default
     * equation and stepper types.
     *  @param[in] fieldType The field type-ID.
     */
    void SetFieldType(G4FieldType fieldType);

    /**
     * Sets or resets the global field. It updates the field objects,
     * if the field was already constructed.
     *  @param[in] field Pointer to the global field.
     *  @param[in] warn If flag is true, issues a warning if the previous
     *             field is deleted.
     */
    void SetGlobalField(G4Field* field, G4bool warn=false);

    /**
     * Registers the local field in the map. It updates the field objects,
     * if the field was already constructed.
     * The field is propagated to all volume daughters regardless if they
     * have already assigned a field manager or not.
     * When multiple local fields are defined (by calling this method multiple
     * times), they will be applied in the order they were set.
     *  @param[in] field Pointer to the global field.
     *  @param[in] lv Pointer to the logical volume.
     *  @param[in] warn If flag is true, issues a warning if the previous
     *             field is deleted.
     */
    void SetLocalField(G4Field* field, G4LogicalVolume* lv, G4bool warn=false);

    /**
     * Sets the user equation of motion.
     *  @param[in] equation Pointer to the equation of motion algorithm.
     *  @param[in] volumeName Optional volume name.
     */
    void SetUserEquationOfMotion(G4EquationOfMotion* equation,
                                 const G4String& volumeName = "");

    /**
     * Sets the user stepper.
     *  @param[in] stepper Pointer to the stepper algorithm.
     *  @param[in] volumeName Optional volume name.
     */
    void SetUserStepper(G4MagIntegratorStepper* stepper,
                        const G4String& volumeName = "");

    /**
     * Sets the verbosity level.
     */
    void SetVerboseLevel(G4int value);

    // Get methods

    /**
     * Gets a pointer to the field parameters with the given 'volumeName'.
     * Return global field parameters, if volume name is empty.
     */
    G4FieldParameters* GetFieldParameters(const G4String& volumeName = "") const;

  private:

    /**
     * Private constructor.
     */
    G4FieldBuilder();

    /**
     * Gets the pointer to field parameters with the given 'volumeName'
     * or creates them if they do not exist yet.
     */
    G4FieldParameters* GetOrCreateFieldParameters(const G4String& volumeName);

    /**
     * Gets the pointer to the field setup with the given logical volume.
     */
    G4FieldSetup* GetFieldSetup(G4LogicalVolume* lv);

    /**
     * Creates magnetic, electromagnetic or gravity field setup.
     */
    void CreateFieldSetup(G4Field* field, G4FieldParameters* fieldParameters,
                          G4LogicalVolume* lv);

    /**
     * Constructs global magnetic field setup.
     */
    void ConstructGlobalField();

    /**
     * Constructs local magnetic field setups from the local fields map.
     */
    void ConstructLocalFields();

    /**
     * Updates all field setups.
     */
    void UpdateFieldSetups();

    /**
     * Helper methods.
     */
    inline std::vector<G4FieldSetup*>& GetFieldSetups();
    inline std::vector<std::pair<G4LogicalVolume*, G4Field*>>& GetLocalFields();

  private:  // Data members

    /** Information if an instance exists. */
    inline static G4ThreadLocal G4bool fgIsInstance { false };

    /** Messenger for this class. */
    G4FieldBuilderMessenger* fMessenger = nullptr;

    /** Field parameters. */
    std::vector<G4FieldParameters*> fFieldParameters;

    /** Field setups. */
    G4Cache<std::vector<G4FieldSetup*>*> fFieldSetups;

    /** Registered global field. */
    static G4ThreadLocal G4Field* fGlobalField;

    /** Registered local fields. */
    G4Cache<std::vector<std::pair<G4LogicalVolume*, G4Field*>>*> fLocalFields;

    /** Info if field objects were constructed. */
    static G4ThreadLocal G4bool fIsConstructed;

    /** Verbose level. */
    G4int fVerboseLevel = 1;
};

// Inline methods -------------------------------------------------------------

inline std::vector<G4FieldSetup*>& G4FieldBuilder::GetFieldSetups()
{
  // Return reference to field setups from G4Cache
  return *fFieldSetups.Get();
}

inline std::vector<std::pair<G4LogicalVolume*, G4Field*>>& G4FieldBuilder::GetLocalFields()
{
  // Return reference to local fields map from G4Cache
  return *fLocalFields.Get();
}

#endif
