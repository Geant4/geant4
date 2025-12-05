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
// G4FieldParametersMessenger
//
// Class description:
//
// Messenger class that defines commands for field configuration.
//
// Implements commands:
// - /field/fieldType fieldType
//       fieldType = Magnetic | ElectroMagnetic | Gravity
// - /field/equationType eqType
//       eqType = EqMagnetic | EqMagneticWithSpin | EqElectroMagnetic |
//                EqEMfieldWithSpin | EqEMfieldWithEDM
// - /field/stepperType stepperType
//       stepperType = CashKarpRKF45 | ClassicalRK4 | ExplicitEuler | ImplicitEuler |
//                     SimpleHeum | SimpleRunge | ConstRK4 | ExactHelixStepper
//                     | HelixExplicitEuler | HelixHeum | HelixImplicitEuler |
//                     HelixMixedStepper | HelixSimpleRunge | NystromRK4 |
//                     RKG3Stepper
// - /field/setMinimumStep value
// - /field/setDeltaChord  value
// - /field/setDeltaOneStep value
// - /field/setDeltaIntersection value
// - /field/setMinimumEpsilonStep value
// - /field/setMaximumEpsilonStep value
// - /field/setConstDistance value
// - /field/printParameters
//
// Only equation type and stepper type values that are handled by G4FieldBuilder
// are accepted by the commands.

// Author: Ivana Hrivnacova (IJClab, Orsay), 2024.
// -------------------------------------------------------------------
#ifndef G4FIELDPARAMETERSMESSENGER_HH
#define G4FIELDPARAMETERSMESSENGER_HH

#include "G4UImessenger.hh"
#include "globals.hh"

class G4FieldParameters;

class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithABool;

/**
 * @brief G4FieldParametersMessenger is a messenger class that defines
 * commands for field configuration. Only equation type and stepper type
 * values that are handled by G4FieldBuilder are accepted by the commands.
 */

class G4FieldParametersMessenger : public G4UImessenger
{
  public:

    /**
     * Standard constructor for G4FieldParametersMessenger.
     *  @param[in] fieldParameters Pointer to the field parameters object.
     */
    G4FieldParametersMessenger(G4FieldParameters* fieldParameters);

    /**
     * Destructor.
     */
    ~G4FieldParametersMessenger() override;

    /**
     * Default constructor, copy constructor and assignment operator not allowed.
     */
    G4FieldParametersMessenger() = delete;
    G4FieldParametersMessenger(const G4FieldParametersMessenger&) = delete;
    G4FieldParametersMessenger& operator=(const G4FieldParametersMessenger&) = delete;

    /**
     * Applies command to the associated object.
     */
    void SetNewValue(G4UIcommand* command, G4String newValues) override;

  private:

    /** Associated class object. */
    G4FieldParameters* fFieldParameters = nullptr;

    /** Commands directory. */
    G4UIdirectory* fDirectory = nullptr;

    // Commands data members

    /** Command: fieldType. */
    G4UIcmdWithAString* fFieldTypeCmd = nullptr; 

    /** Command: equationType. */
    G4UIcmdWithAString* fEquationTypeCmd = nullptr; 

    /** Command: stepperType. */
    G4UIcmdWithAString* fStepperTypeCmd = nullptr; 

    /** Command: setMinimumStep. */
    G4UIcmdWithADoubleAndUnit* fSetMinimumStepCmd = nullptr; 

    /** Command: setDeltaChord. */
    G4UIcmdWithADoubleAndUnit* fSetDeltaChordCmd = nullptr; 

    /** Command: setDeltaOneStep. */
    G4UIcmdWithADoubleAndUnit* fSetDeltaOneStepCmd = nullptr; 

    /** Command: setDeltaIntersection. */
    G4UIcmdWithADoubleAndUnit* fSetDeltaIntersectionCmd = nullptr; 

    /** Command: setMinimumEpsilon. */
    G4UIcmdWithADouble* fSetMinimumEpsilonStepCmd = nullptr; 

    /** Command: setMaximumEpsilon. */
    G4UIcmdWithADouble* fSetMaximumEpsilonStepCmd = nullptr; 

    /** Command: setConstDistance. */
    G4UIcmdWithADoubleAndUnit* fSetConstDistanceCmd = nullptr; 

    /** Command: printParameters. */
    G4UIcmdWithoutParameter* fPrintParametersCmd = nullptr; 
};

#endif
