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
// G4FieldParameters
//
// Class description:
//
// The class defines the type of equation of motion of a particle
// in a field and the integration method, as well as other accuracy
// parameters.
//
// The default values correspond to the defaults set in Geant4.

// Author: Ivana Hrivnacova (IJCLab, Orsay), 2024.
// -------------------------------------------------------------------
#ifndef G4FIELDPARAMETERS_HH
#define G4FIELDPARAMETERS_HH

#include "globals.hh"

#include <CLHEP/Units/SystemOfUnits.h>

class G4FieldParametersMessenger;
class G4EquationOfMotion;
class G4MagIntegratorStepper;

/**
 * @brief G4FieldType defines the available fields in Geant4.
 */

enum G4FieldType
{
  kMagnetic,           ///< magnetic field
  kElectroMagnetic,    ///< electromagnetic field
  kGravity,            ///< gravity field
  kUserFieldType       ///< User defined field type
};

/**
 * @brief G4EquationType defines the types of equations of motion of a
 * particle in a field in Geant4.
 */

enum G4EquationType
{
  kEqMagnetic,        ///< G4Mag_UsualEqRhs: the standard right-hand side for
                      ///< equation of motion.
  kEqMagneticWithSpin,///< G4Mag_SpinEqRhs: the equation of motion for a particle
                      ///< with spin
                      ///< in a pure magnetic field
  kEqElectroMagnetic, ///< G4EqMagElectricField: Equation of motion in a combined
                      ///< electric and magnetic field
  kEqEMfieldWithSpin, ///< G4EqEMFieldWithSpin: Equation of motion for a
                      ///< particle with spin
                      ///< in a combined electric and magnetic field
  kEqEMfieldWithEDM,  ///< G4EqEMFieldWithEDM: Equation of motion in a combined
                      ///< electric and magnetic field, with spin tracking for
                      ///< both MDM and EDM terms
  kEqGravity,         ///< G4EqGravityField: equation of motion in a gravity field
                      ///  (not build by G4FieldBuilder)
  kEqMonopole,        ///< G4MonopoleEq: the right-hand side of equation of motion for monopole
                      ///  in a combined electric and magnetic field
                      ///  (not build by G4FieldBuilder)
  kEqReplate,         ///< G4RepleteEofM: equation of motion in a combined field, including:
                      ///  magnetic, electric, gravity, and gradient B field, as well as spin tracking
                      ///  (not build by G4FieldBuilder)
  kUserEquation       ///< User defined equation of motion
};

/**
 * @brief G4StepperType defines the available integrator of particle's
 * equation of motion in Geant4.
 */

enum G4StepperType
{
  // steppers with equation of motion of generic type (G4EquationOfMotion)
  kCashKarpRKF45,     ///< G4CashKarpRKF45
  kClassicalRK4,      ///< G4ClassicalRK4
  kBogackiShampine23, ///< G4BogackiShampine23
  kBogackiShampine45, ///< G4BogackiShampine45
  kDoLoMcPriRK34,     ///< G4DoLoMcPriRK34
  kDormandPrince745,  ///< G4DormandPrince745
  kDormandPrinceRK56, ///< G4DormandPrinceRK56
  kDormandPrinceRK78, ///< G4DormandPrinceRK78
  kExplicitEuler,     ///< G4ExplicitEuler
  kImplicitEuler,     ///< G4ImplicitEuler
  kSimpleHeum,        ///< G4SimpleHeum
  kSimpleRunge,       ///< G4SimpleRunge
  kTsitourasRK45,     ///< G4TsitourasRK45

  // steppers with equation of motion of G4Mag_UsualEqRhs type
  kConstRK4,           ///< G4ConstRK4
  kExactHelixStepper,  ///< G4ExactHelixStepper
  kHelixExplicitEuler, ///< G4HelixExplicitEuler
  kHelixHeum,          ///< G4HelixHeum
  kHelixImplicitEuler, ///< G4HelixImplicitEuler
  kHelixMixedStepper,  ///< G4HelixMixedStepper
  kHelixSimpleRunge,   ///< G4HelixSimpleRunge
  kNystromRK4,         ///< G4NystromRK4
  kRKG3Stepper,        ///< G4RKG3_Stepper
  kUserStepper,        ///< User defined stepper

  // FSAL steppers
  kRK547FEq1,          ///< G4RK547FEq1
  kRK547FEq2,          ///< G4RK547FEq2
  kRK547FEq3,          ///< G4RK547FEq3

  // Templated steppers (not build by G4FieldBuilder)
  kTCashKarpRKF45,     ///< G4TCashKarpRKF45
  kTDormandPrince45,   ///< G4TDormandPrince45
  kTMagErrorStepper,   ///< G4TMagErrorStepper
  kQSStepper           ///< G4QSStepper
};

/**
 * @brief G4FieldDefaults defines the magnetic field parameters defaults.
 * The namespace defines the default values of the field paraments as constexpr
 * so that they can be used also as the default values in the magnetic field
 * classes constructors and other member functions.
 */

namespace G4FieldDefaults
{
  /// Default minimum step in G4ChordFinder
  constexpr G4double kMinimumStep  = 0.01 * CLHEP::mm;
  /// Default delta chord in G4ChordFinder
  constexpr G4double kDeltaChord = 0.25 * CLHEP::mm;
  /// Default delta one step in global field manager
  constexpr G4double kDeltaOneStep = 0.01 * CLHEP::mm;
  /// Delta intersection in global field manager
  constexpr G4double kDeltaIntersection = 0.001 * CLHEP::mm;
  /// Default minimum epsilon step in global field manager
  constexpr G4double kMinimumEpsilonStep = 5.0e-5;  // Expected: 5.0e-5 to 1.0e-10 ...
  /// Default maximum epsilon step in global field manager
  constexpr G4double kMaximumEpsilonStep = 0.001;   // Expected: 1.0e-3 to 1.0e-8 ...
}

/**
 * @brief G4FieldParameters defines the type of equation of motion of a
 * particle in a field and the integration method, as well as other accuracy
 * parameters. The default values correspond to the defaults set in Geant4.
 */

class G4FieldParameters
{
  public:

    /**
     * Constructor for G4FieldParameters.
     *  @param[in] volumeName The volume name where field is applied.
     */
    G4FieldParameters(const G4String& volumeName = "");

    /**
     * Destructor.
     */
    ~G4FieldParameters();

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4FieldParameters(const G4FieldParameters& right) = delete;
    G4FieldParameters& operator=(const G4FieldParameters& right) = delete;

    /**
     * Returns the field type as a string.
     */
    static G4String FieldTypeName(G4FieldType field);

    /**
     * Returns the equation type as a string.
     */
    static G4String EquationTypeName(G4EquationType equation);

    /**
     * Returns the stepper type as a string.
     */
    static G4String StepperTypeName(G4StepperType stepper);

    /**
     * Returns the field type for given field type name.
     */
    static G4FieldType GetFieldType(const G4String& name);

    /**
     * Returns the equation type for given equation type name.
     */
    static G4EquationType GetEquationType(const G4String& name);

    /**
     * Returns the stepper type for given stepper type name.
     */
    static G4StepperType GetStepperType(const G4String& name);

    /**
     * Prints all customisable accuracy parameters.
     */
    void PrintParameters() const;

    // Set methods ------------------------------------------------------------

    /**
     * Sets the type of field.
     */
    void SetFieldType(G4FieldType field);

    /**
     * Sets the type of equation of motion of a particle in a field.
     */
    void SetEquationType(G4EquationType equation);

    /**
     * Sets the type of integrator of particle's equation of motion.
     */
    void SetStepperType(G4StepperType stepper);

    /**
     * Sets the user defined equation of motion.
     */
    void SetUserEquationOfMotion(G4EquationOfMotion* equation);

    /**
     * Sets the user defined integrator of particle's equation of motion.
     */
    void SetUserStepper(G4MagIntegratorStepper* stepper);

    /**
     * Sets the minimum step in G4ChordFinder.
     */
    void SetMinimumStep(G4double value);

    /**
     * Sets the delta chord in G4ChordFinder.
     */
    void SetDeltaChord(G4double value);

    /**
     * Sets the delta one step in global field manager.
     */
    void SetDeltaOneStep(G4double value);

    /**
     * Sets the delta intersection in global field manager.
     */
    void SetDeltaIntersection(G4double value);

    /**
     * Sets the minimum epsilon step in global field manager.
     */
    void SetMinimumEpsilonStep(G4double value);

    /**
     * Sets the maximum epsilon step in global field manager.
     */
    void SetMaximumEpsilonStep(G4double value);

    /**
     * Sets the distance within which the field is considered constant.
     */
    void SetConstDistance(G4double value);

    // Get methods ------------------------------------------------------------

    /**
     * Gets the name of associated volume, if local field.
     */
    const G4String& GetVolumeName() const;

    /**
     * Gets the type of field.
     */
    const G4FieldType& GetFieldType() const;

    /**
     * Gets the type of equation of motion of a particle in a field.
     */
    const G4EquationType& GetEquationType() const;

    /**
     * Gets the type of integrator of particle's equation of motion.
     */
    const G4StepperType& GetStepperType() const;

    /**
     * Gets the user defined equation of motion.
     */
    G4EquationOfMotion* GetUserEquationOfMotion() const;

    /**
     * Gets the user defined integrator of particle's equation of motion.
     */
    G4MagIntegratorStepper* GetUserStepper() const;

    /**
     * Gets the minimum step in G4ChordFinder.
     */
    G4double GetMinimumStep() const;

    /**
     * Gets the delta chord in G4ChordFinder.
     */
    G4double GetDeltaChord() const;

    /**
     * Gets the delta one step in global field manager.
     */
    G4double GetDeltaOneStep() const;

    /**
     * Gets the delta intersection in global field manager.
     */
    G4double GetDeltaIntersection() const;

    /**
     * Gets the minimum epsilon step in global field manager.
     */
    G4double GetMinimumEpsilonStep() const;

    /**
     * Gets the maximum epsilon step in global field manager.
     */
    G4double GetMaximumEpsilonStep() const;

    /**
     * Gets the distance within which the field is considered constant.
     */
    G4double GetConstDistance() const;

  private:

    /** Default constant distance. */
    inline static const G4double fgkDefaultConstDistance = 0.;

    /** Messenger for this class. */
    G4FieldParametersMessenger* fMessenger = nullptr;

    /** The name of the associated volume, if local field. */
    G4String fVolumeName;

    /** The minimum step in G4ChordFinder. */
    G4double fMinimumStep = G4FieldDefaults::kMinimumStep;

    /** The delta chord in G4ChordFinder. */
    G4double fDeltaChord = G4FieldDefaults::kDeltaChord;

    /** The delta one step in global field manager. */
    G4double fDeltaOneStep = G4FieldDefaults::kDeltaOneStep;

    /** The delta intersection in global field manager. */
    G4double fDeltaIntersection = G4FieldDefaults::kDeltaIntersection;

    /** The minimum epsilon step in global field manager. */
    G4double fMinimumEpsilonStep = G4FieldDefaults::kMinimumEpsilonStep;

    /** The maximum epsilon step in global field manager. */
    G4double fMaximumEpsilonStep = G4FieldDefaults::kMaximumEpsilonStep;

    /** The type of field. */
    G4FieldType fField = kMagnetic;

    /** Type of equation of motion of a particle in a field. */
    G4EquationType fEquation = kEqMagnetic;

    /** Type of integrator of particle's equation of motion. */
    G4StepperType fStepper = kDormandPrince745;

    /** User defined equation of motion. */
    G4EquationOfMotion* fUserEquation = nullptr;

    /// User defined integrator of particle's equation of motion. */
    G4MagIntegratorStepper* fUserStepper = nullptr;

    /** The distance within which the field is considered constant. */
    G4double fConstDistance = fgkDefaultConstDistance;
};

// Inline functions

#include "G4FieldParameters.icc"

#endif
