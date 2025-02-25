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

/// \file G4FieldParameters.cc
/// \brief Implementation of the G4FieldParameters class
///
/// This code was initially developed in Geant4 VMC package
/// (https://github.com/vmc-project)
/// and adapted to Geant4.
///
/// \author I. Hrivnacova; IJCLab, Orsay

#include "G4FieldParameters.hh"
#include "G4FieldParametersMessenger.hh"

#include "G4Exception.hh"
#include "G4SystemOfUnits.hh"

//
// static methods
//

//_____________________________________________________________________________
G4String G4FieldParameters::FieldTypeName(G4FieldType field)
{
  // Return the field type as a string

  switch (field) {
    case kMagnetic:
      return G4String("Magnetic");
    case kElectroMagnetic:
      return G4String("ElectroMagnetic");
    case kGravity:
      return G4String("Gravity");
  }

  G4Exception(
    "G4FieldParameters::FieldTypeName:", "GeomFieldParameters0001",
    FatalErrorInArgument, "Unknown field value.");
  return G4String();
}

//_____________________________________________________________________________
G4String G4FieldParameters::EquationTypeName(G4EquationType equation)
{
  // Return the equation type as a string

  switch (equation) {
    case kEqMagnetic:
      return G4String("EqMagnetic");
    case kEqMagneticWithSpin:
      return G4String("EqMagneticWithSpin");
    case kEqElectroMagnetic:
      return G4String("EqElectroMagnetic");
    case kEqEMfieldWithSpin:
      return G4String("EqEMfieldWithSpin");
    case kEqEMfieldWithEDM:
      return G4String("EqEMfieldWithEDM");
    case kUserEquation:
      return G4String("UserDefinedEq");
  }

  G4Exception(
    "G4FieldParameters::EquationTypeName:", "GeomFieldParameters0001",
    FatalErrorInArgument, "Unknown equation value.");
  return G4String();
}

//_____________________________________________________________________________
G4String G4FieldParameters::StepperTypeName(G4StepperType stepper)
{
  // Return the stepper type as a string

  switch (stepper) {
    case kBogackiShampine23:
      return G4String("BogackiShampine23");
    case kBogackiShampine45:
      return G4String("BogackiShampine45");
    case kCashKarpRKF45:
      return G4String("CashKarpRKF45");
    case kClassicalRK4:
      return G4String("ClassicalRK4");
    case kDormandPrince745:
      return G4String("DormandPrince745");
    case kDormandPrinceRK56:
      return G4String("DormandPrinceRK56");
    case kDormandPrinceRK78:
      return G4String("DormandPrinceRK78");
    case kExplicitEuler:
      return G4String("ExplicitEuler");
    case kImplicitEuler:
      return G4String("ImplicitEuler");
    case kSimpleHeum:
      return G4String("SimpleHeum");
    case kSimpleRunge:
      return G4String("SimpleRunge");
    case kConstRK4:
      return G4String("ConstRK4");
    case kExactHelixStepper:
      return G4String("ExactHelixStepper");
    case kHelixExplicitEuler:
      return G4String("HelixExplicitEuler");
    case kHelixHeum:
      return G4String("HelixHeum");
    case kHelixImplicitEuler:
      return G4String("HelixImplicitEuler");
    case kHelixMixedStepper:
      return G4String("HelixMixedStepper");
    case kHelixSimpleRunge:
      return G4String("HelixSimpleRunge");
    case kNystromRK4:
      return G4String("NystromRK4");
    case kRKG3Stepper:
      return G4String("RKG3_Stepper");
    case kTsitourasRK45:
      return G4String("TsitourasRK45");
    case kUserStepper:
      return G4String("UserDefinedStepper");
    case kRK547FEq1:
      return G4String("RK547FEq1");
    case kRK547FEq2:
      return G4String("RK547FEq2");
    case kRK547FEq3:
      return G4String("RK547FEq3");
  }

  G4Exception(
    "G4FieldParameters::StepperTypeName:", "GeomFieldParameters0001",
    FatalErrorInArgument, "Unknown stepper value.");
  return G4String();
}

//_____________________________________________________________________________
G4FieldType G4FieldParameters::GetFieldType(const G4String& name)
{
  // Return the field type for given field type name

  if (name == FieldTypeName(kMagnetic)) return kMagnetic;
  if (name == FieldTypeName(kElectroMagnetic)) return kElectroMagnetic;
  if (name == FieldTypeName(kGravity)) return kGravity;

  G4Exception(
    "G4FieldParameters::GetFieldType:", "GeomFieldParameters0001",
    FatalErrorInArgument, "Unknown field name.");
  return kMagnetic;
}

//_____________________________________________________________________________
G4EquationType G4FieldParameters::GetEquationType(const G4String& name)
{
  // Return the equation type for given equation type name

  if (name == EquationTypeName(kEqMagnetic)) return kEqMagnetic;
  if (name == EquationTypeName(kEqMagneticWithSpin)) return kEqMagneticWithSpin;
  if (name == EquationTypeName(kEqElectroMagnetic)) return kEqElectroMagnetic;
  if (name == EquationTypeName(kEqEMfieldWithSpin)) return kEqEMfieldWithSpin;
  if (name == EquationTypeName(kEqEMfieldWithEDM)) return kEqEMfieldWithEDM;
  if (name == EquationTypeName(kUserEquation)) return kUserEquation;

  G4Exception(
    "G4FieldParameters::GetEquationType:", "GeomFieldParameters0001",
    FatalErrorInArgument, "Unknown equation name.");
  return kEqMagnetic;
}

//_____________________________________________________________________________
G4StepperType G4FieldParameters::GetStepperType(const G4String& name)
{
  // Return the stepper type for given stepper type name
  if (name == StepperTypeName(kBogackiShampine23)) return kBogackiShampine23;
  if (name == StepperTypeName(kBogackiShampine45)) return kBogackiShampine45;
  if (name == StepperTypeName(kCashKarpRKF45)) return kCashKarpRKF45;
  if (name == StepperTypeName(kClassicalRK4)) return kClassicalRK4;
  if (name == StepperTypeName(kDormandPrince745)) return kDormandPrince745;
  if (name == StepperTypeName(kDormandPrinceRK56)) return kDormandPrinceRK56;
  if (name == StepperTypeName(kDormandPrinceRK78)) return kDormandPrinceRK78;
  if (name == StepperTypeName(kExplicitEuler)) return kExplicitEuler;
  if (name == StepperTypeName(kImplicitEuler)) return kImplicitEuler;
  if (name == StepperTypeName(kSimpleHeum)) return kSimpleHeum;
  if (name == StepperTypeName(kSimpleRunge)) return kSimpleRunge;
  if (name == StepperTypeName(kConstRK4)) return kConstRK4;
  if (name == StepperTypeName(kExactHelixStepper)) return kExactHelixStepper;
  if (name == StepperTypeName(kHelixExplicitEuler)) return kHelixExplicitEuler;
  if (name == StepperTypeName(kHelixHeum)) return kHelixHeum;
  if (name == StepperTypeName(kHelixImplicitEuler)) return kHelixImplicitEuler;
  if (name == StepperTypeName(kHelixMixedStepper)) return kHelixMixedStepper;
  if (name == StepperTypeName(kHelixSimpleRunge)) return kHelixSimpleRunge;
  if (name == StepperTypeName(kNystromRK4)) return kNystromRK4;
  if (name == StepperTypeName(kRKG3Stepper)) return kRKG3Stepper;
  if (name == StepperTypeName(kRK547FEq1)) return kRK547FEq1;
  if (name == StepperTypeName(kRK547FEq2)) return kRK547FEq2;
  if (name == StepperTypeName(kRK547FEq3)) return kRK547FEq3;
  if (name == StepperTypeName(kTsitourasRK45)) return kTsitourasRK45;
  if (name == StepperTypeName(kUserStepper)) return kUserStepper;

  G4Exception(
    "G4FieldParameters::GetStepperType:", "GeomFieldParameters0001",
    FatalErrorInArgument, "Unknown stepper name.");
  return kClassicalRK4;
}

//
// ctors, dtor
//

//_____________________________________________________________________________
G4FieldParameters::G4FieldParameters(const G4String& volumeName)
  : fVolumeName(volumeName)
{
  // Standard constructor

  fMessenger = new G4FieldParametersMessenger(this);
}

//_____________________________________________________________________________
G4FieldParameters::~G4FieldParameters()
{
  // Destructor

  delete fMessenger;
}

//
// public methods
//

//_____________________________________________________________________________
void G4FieldParameters::PrintParameters() const
{
  // Prints all customizable accuracy parameters

  G4cout << "Magnetic field parameters: " << G4endl;
  if (fVolumeName.size()) {
    G4cout << "  volume name = " << fVolumeName << G4endl;
  }
  G4cout << "  field type = " << FieldTypeName(fField) << G4endl
         << "  equation type = " << EquationTypeName(fEquation) << G4endl
         << "  stepper type = " << StepperTypeName(fStepper) << G4endl
         << "  minStep = " << fMinimumStep << " mm" << G4endl
         << "  constDistance = " << fConstDistance << " mm" << G4endl
         << "  deltaChord = " << fDeltaChord << " mm" << G4endl
         << "  deltaOneStep = " << fDeltaOneStep << " mm" << G4endl
         << "  deltaIntersection = " << fDeltaIntersection << " mm" << G4endl
         << "  epsMin = " << fMinimumEpsilonStep << G4endl
         << "  epsMax=  " << fMaximumEpsilonStep << G4endl;
}

//_____________________________________________________________________________
void G4FieldParameters::SetUserEquationOfMotion(G4EquationOfMotion* equation)
{
  // Set user defined equation of motion

  fUserEquation = equation;
  fEquation = kUserEquation;
}

//_____________________________________________________________________________
void G4FieldParameters::SetUserStepper(G4MagIntegratorStepper* stepper)
{
  // Set user defined integrator of particle's equation of motion

  fUserStepper = stepper;
  fStepper = kUserStepper;
}
