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

/// \file G4FieldParameters.hh
/// \brief Definition of the G4FieldParameters class
///
/// This code was initially developed in Geant4 VMC package
/// (https://github.com/vmc-project)
/// and adapted to Geant4.
///
/// \author I. Hrivnacova; IJCLab, Orsay

#ifndef G4FIELDPARAMETERS_HH
#define G4FIELDPARAMETERS_HH

#include "G4MagneticField.hh"
#include "globals.hh"

#include <CLHEP/Units/SystemOfUnits.h>

class G4FieldParametersMessenger;

class G4EquationOfMotion;
class G4MagIntegratorStepper;

/// The available fields in Geant4
enum G4FieldType
{
  kMagnetic,        ///< magnetic field
  kElectroMagnetic, ///< electromagnetic field
  kGravity          ///< gravity field
};

/// The available equations of motion of a particle in a field
/// in Geant4
enum G4EquationType
{
  kMagUsualEqRhs,     ///< G4Mag_UsualEqRhs: the standard right-hand side for
                      ///< equation
                      /// of motion.
  kMagSpinEqRhs,      ///< G4Mag_SpinEqRhs: the equation of motion for a particle
                      ///< with spin
                      /// in a pure magnetic field
  kEqMagElectric,     ///< G4EqMagElectricField: Equation of motion in a combined
                      ///  electric and magnetic field
  kEqEMFieldWithSpin, ///< G4EqEMFieldWithSpin: Equation of motion for a
                      ///< particle with spin
                      ///  in a combined electric and magnetic field
  kEqEMFieldWithEDM,  ///< G4EqEMFieldWithEDM: Equation of motion in a combined
                      /// electric and magnetic field, with spin tracking for
                      /// both MDM and EDM terms
  kUserEquation       ///< User defined equation of motion
};

/// The available integrator of particle's equation of motion
/// in Geant4
enum G4StepperType
{
  // steppers with equation of motion of generic type (G4EquationOfMotion)
  kCashKarpRKF45,     ///< G4CashKarpRKF45
  kClassicalRK4,      ///< G4ClassicalRK4
  kBogackiShampine23, ///< G4BogackiShampine23
  kBogackiShampine45, ///< G4BogackiShampine45
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
  kRK547FEq1, ///< G4RK547FEq1
  kRK547FEq2, ///< G4RK547FEq2
  kRK547FEq3  ///< G4RK547FEq3
};

/// \brief The magnetic field parameters
///
/// The class defines the type of equation of motion of a particle
/// in a field and the integration method, as well as other accuracy
/// parameters.
///
/// The default values correspond to the defaults set in Geant4
/// (taken from Geant4 9.3 release.)
/// As Geant4 classes to not provide access methods for these defaults,
/// the defaults have to be checked with each new Geant4 release.
///
/// \author I. Hrivnacova; IJCLab, Orsay

class G4FieldParameters
{
 public:
  /// Standard and default constructor
  G4FieldParameters(const G4String& volumeName = "");
  /// Destructor
  ~G4FieldParameters();

  // Methods
  //

  /// Return the field type as a string
  static G4String FieldTypeName(G4FieldType field);
  /// Return the equation type as a string
  static G4String EquationTypeName(G4EquationType equation);
  /// Return the stepper type as a string
  static G4String StepperTypeName(G4StepperType stepper);
  /// Return the field type for given field type name
  static G4FieldType GetFieldType(const G4String& name);
  /// Return the equation type for given equation type name
  static G4EquationType GetEquationType(const G4String& name);
  /// Return the stepper type for given stepper type name
  static G4StepperType GetStepperType(const G4String& name);

  /// Prints all customizable accuracy parameters
  void PrintParameters() const;

  // Set methods
  //

  /// Set type of field
  void SetFieldType(G4FieldType field);
  /// Set Type of equation of motion of a particle in a field
  void SetEquationType(G4EquationType equation);
  /// Type of integrator of particle's equation of motion
  void SetStepperType(G4StepperType stepper);
  /// Set user defined equation of motion
  void SetUserEquationOfMotion(G4EquationOfMotion* equation);
  /// Set user defined integrator of particle's equation of motion
  void SetUserStepper(G4MagIntegratorStepper* stepper);

  /// Set minimum step in G4ChordFinder
  void SetMinimumStep(G4double value);
  /// Set delta chord in G4ChordFinder
  void SetDeltaChord(G4double value);
  /// Set delta one step in global field manager
  void SetDeltaOneStep(G4double value);
  /// Set delta intersection in global field manager
  void SetDeltaIntersection(G4double value);
  /// Set minimum epsilon step in global field manager
  void SetMinimumEpsilonStep(G4double value);
  /// Set maximum epsilon step in global field manager
  void SetMaximumEpsilonStep(G4double value);
  /// Set the distance within which the field is considered constant
  void SetConstDistance(G4double value);

  // Get methods
  //

  // Get the name of associated volume, if local field
  G4String GetVolumeName() const;

  /// Get type of field
  G4FieldType GetFieldType() const;
  /// Get type of equation of motion of a particle in a field
  G4EquationType GetEquationType() const;
  /// Get rype of integrator of particle's equation of motion
  G4StepperType GetStepperType() const;
  /// Get user defined equation of motion
  G4EquationOfMotion* GetUserEquationOfMotion() const;
  /// Get user defined integrator of particle's equation of motion
  G4MagIntegratorStepper* GetUserStepper() const;

  /// Get minimum step in G4ChordFinder
  G4double GetMinimumStep() const;
  /// Get delta chord in G4ChordFinder
  G4double GetDeltaChord() const;
  /// Get delta one step in global field manager
  G4double GetDeltaOneStep() const;
  /// Get delta intersection in global field manager
  G4double GetDeltaIntersection() const;
  /// Get minimum epsilon step in global field manager
  G4double GetMinimumEpsilonStep() const;
  /// Get maximum epsilon step in global field manager
  G4double GetMaximumEpsilonStep() const;
  /// Get the distance within which the field is considered constant
  G4double GetConstDistance() const;

 private:
  /// Not implemented
  G4FieldParameters(const G4FieldParameters& right) = delete;
  /// Not implemented
  G4FieldParameters& operator=(const G4FieldParameters& right) = delete;

  // static data members
  //
  /// Default minimum step in G4ChordFinder
  inline static const G4double fgkDefaultMinimumStep  = 0.01 * CLHEP::mm;
  /// Default delta chord in G4ChordFinder
  inline static const G4double fgkDefaultDeltaChord = 0.25 * CLHEP::mm;
  /// Default delta one step in global field manager
  inline static const G4double fgkDefaultDeltaOneStep = 0.01 * CLHEP::mm;
  /// Delta intersection in global field manager
  inline static const G4double fgkDefaultDeltaIntersection = 0.001 * CLHEP::mm;
  /// Default minimum epsilon step in global field manager
  inline static const G4double fgkDefaultMinimumEpsilonStep = 5.0e-5;
  /// Default maximum epsilon step in global field manager
  inline static const G4double fgkDefaultMaximumEpsilonStep = 0.001;
  /// Default constant distance
  inline static const G4double fgkDefaultConstDistance = 0.;

  // data members
  //
  /// Messenger for this class
  G4FieldParametersMessenger* fMessenger = nullptr;

  /// The name of associated volume, if local field
  G4String fVolumeName;

  /// Minimum step in G4ChordFinder
  G4double fMinimumStep = fgkDefaultMinimumStep;
  /// Delta chord in G4ChordFinder
  G4double fDeltaChord = fgkDefaultDeltaChord;
  /// Delta one step in global field manager
  G4double fDeltaOneStep = fgkDefaultDeltaOneStep;
  /// Delta intersection in global field manager
  G4double fDeltaIntersection = fgkDefaultDeltaIntersection;
  /// Minimum epsilon step in global field manager
  G4double fMinimumEpsilonStep = fgkDefaultMinimumEpsilonStep;
  /// Maximum epsilon step in global field manager
  G4double fMaximumEpsilonStep = fgkDefaultMaximumEpsilonStep;

  /// Type of field
  G4FieldType fField = kMagnetic;

  /// Type of equation of motion of a particle in a field
  G4EquationType fEquation = kMagUsualEqRhs;

  /// Type of integrator of particle's equation of motion
  G4StepperType fStepper = kDormandPrince745;

  /// User defined equation of motion
  G4EquationOfMotion* fUserEquation = nullptr;

  /// User defined integrator of particle's equation of motion
  G4MagIntegratorStepper* fUserStepper = nullptr;

  /// The distance within which the field is considered constant
  G4double fConstDistance = fgkDefaultConstDistance;
};

// inline functions

// Set type of field
inline void G4FieldParameters::SetFieldType(G4FieldType field)
{
  fField = field;
}

// Set the type of equation of motion of a particle in a field
inline void G4FieldParameters::SetEquationType(G4EquationType equation)
{
  fEquation = equation;
}

// Set the type of integrator of particle's equation of motion
inline void G4FieldParameters::SetStepperType(G4StepperType stepper)
{
  fStepper = stepper;
}

// Set minimum step in G4ChordFinder
inline void G4FieldParameters::SetMinimumStep(G4double value)
{
  fMinimumStep = value;
}

// Set delta chord in G4ChordFinder
inline void G4FieldParameters::SetDeltaChord(G4double value)
{
  fDeltaChord = value;
}

// Set delta one step in global field manager
inline void G4FieldParameters::SetDeltaOneStep(G4double value)
{
  fDeltaOneStep = value;
}

// Set delta intersection in global field manager
inline void G4FieldParameters::SetDeltaIntersection(G4double value)
{
  fDeltaIntersection = value;
}

// Set minimum epsilon step in global field manager
inline void G4FieldParameters::SetMinimumEpsilonStep(G4double value)
{
  fMinimumEpsilonStep = value;
}

// Set maximum epsilon step in global field manager
inline void G4FieldParameters::SetMaximumEpsilonStep(G4double value)
{
  fMaximumEpsilonStep = value;
}

// Set the distance within which the field is considered constant
inline void G4FieldParameters::SetConstDistance(G4double value)
{
  fConstDistance = value;
}

// Return the name of associated volume, if local field
inline G4String G4FieldParameters::GetVolumeName() const
{
  return fVolumeName;
}

// Return the type of field
inline G4FieldType G4FieldParameters::GetFieldType() const { return fField; }

// Return the type of equation of motion of a particle in a field
inline G4EquationType G4FieldParameters::GetEquationType() const
{
  return fEquation;
}

// Return the type of integrator of particle's equation of motion
inline G4StepperType G4FieldParameters::GetStepperType() const
{
  return fStepper;
}

// Return the user defined equation of motion
inline G4EquationOfMotion* G4FieldParameters::GetUserEquationOfMotion() const
{
  return fUserEquation;
}

// Return the user defined integrator of particle's equation of motion
inline G4MagIntegratorStepper* G4FieldParameters::GetUserStepper() const
{
  return fUserStepper;
}

// Return minimum step in G4ChordFinder
inline G4double G4FieldParameters::GetMinimumStep() const
{
  return fMinimumStep;
}

// Return delta chord in G4ChordFinder
inline G4double G4FieldParameters::GetDeltaChord() const
{
  return fDeltaChord;
}

// Return delta one step in global field manager
inline G4double G4FieldParameters::GetDeltaOneStep() const
{
  return fDeltaOneStep;
}

// Return delta intersection in global field manager
inline G4double G4FieldParameters::GetDeltaIntersection() const
{
  return fDeltaIntersection;
}

// Return minimum epsilon step in global field manager
inline G4double G4FieldParameters::GetMinimumEpsilonStep() const
{
  return fMinimumEpsilonStep;
}

// Return maximum epsilon step in global field manager
inline G4double G4FieldParameters::GetMaximumEpsilonStep() const
{
  return fMaximumEpsilonStep;
}

// Return the distance within which the field is considered constant
inline G4double G4FieldParameters::GetConstDistance() const
{
  return fConstDistance;
}

#endif // G4FIELDPARAMETERS_HH
