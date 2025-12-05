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
// Implementation of the G4FieldSetup class
//
// Author: Ivana Hrivnacova (IJClab, Orsay), 2024.
// --------------------------------------------------------------------

#include "G4FieldSetup.hh"
#include "G4FieldSetupMessenger.hh"

#include "G4Exception.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"

#include "G4BogackiShampine23.hh"
#include "G4BogackiShampine45.hh"
#include "G4CachedMagneticField.hh"
#include "G4CashKarpRKF45.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4ConstRK4.hh"
#include "G4DoLoMcPriRK34.hh"
#include "G4DormandPrince745.hh"
#include "G4DormandPrinceRK56.hh"
#include "G4DormandPrinceRK78.hh"
#include "G4ElectroMagneticField.hh"
#include "G4EqEMFieldWithEDM.hh"
#include "G4EqEMFieldWithSpin.hh"
#include "G4EqMagElectricField.hh"
#include "G4ExactHelixStepper.hh"
#include "G4ExplicitEuler.hh"
#include "G4FSALIntegrationDriver.hh"
#include "G4FieldManager.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixHeum.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixMixedStepper.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4ImplicitEuler.hh"
#include "G4MagneticField.hh"
#include "G4MagErrorStepper.hh"
#include "G4MagHelicalStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4Mag_EqRhs.hh"
#include "G4Mag_SpinEqRhs.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4NystromRK4.hh"
#include "G4RK547FEq1.hh"
#include "G4RK547FEq2.hh"
#include "G4RK547FEq3.hh"
#include "G4RKG3_Stepper.hh"
#include "G4SimpleHeum.hh"
#include "G4SimpleRunge.hh"
#include "G4TsitourasRK45.hh"
#include "G4VIntegrationDriver.hh"


//_____________________________________________________________________________
G4FieldSetup::G4FieldSetup(const G4FieldParameters& parameters,
                           G4Field* field, G4LogicalVolume* lv)
  : fParameters(parameters), fG4Field(field), fLogicalVolume(lv)
{
  // Standard constructor

  fMessenger = new G4FieldSetupMessenger(this);

  // Get or create field manager
  if (fLogicalVolume == nullptr)
  {
    // global field
    fFieldManager = G4FieldManager::GetGlobalFieldManager();
  }
  else
  {
    // local field
    fFieldManager = new G4FieldManager();
    G4bool overwriteDaughtersField = true;
      // TO DO: this parameter should be made optional for users
    fLogicalVolume->SetFieldManager(fFieldManager, overwriteDaughtersField);
  }
}

//_____________________________________________________________________________
G4FieldSetup::~G4FieldSetup()
{
  // Destructor
  delete fG4Field;
  delete fChordFinder;
  delete fStepper;
}

//
// private methods
//

//_____________________________________________________________________________
G4Field* G4FieldSetup::CreateCachedField(
    const G4FieldParameters& parameters, G4Field* field)
{
  // Create cached magnetic field if const distance is set > 0.
  // and field is of G4MagneticField.
  // Return the input field otherwise.

  if (parameters.GetConstDistance() > 0.)
  {
    auto magField = dynamic_cast<G4MagneticField*>(field);
    if (magField == nullptr)
    {
      G4Exception(
        "G4FieldSetup::CreateCachedField:", "GeomFieldParameters0001",
        JustWarning, "Incompatible field type.");
      return field;
    }
    return new G4CachedMagneticField(magField, parameters.GetConstDistance());
  }

  return field;
}

//_____________________________________________________________________________
G4EquationOfMotion* G4FieldSetup::CreateEquation(G4EquationType equation)
{
  // Set the equation of motion of a particle in a field

  // magnetic fields
  G4MagneticField* magField = nullptr;
  if (equation == kEqMagnetic || equation == kEqMagneticWithSpin)
  {
    magField = dynamic_cast<G4MagneticField*>(fG4Field);
    if (magField == nullptr)
    {
      G4Exception(
        "G4FieldSetup::CreateEquation:", "GeomFieldParameters0001",
        FatalErrorInArgument, "Incompatible field and equation.\n"
        "The field type must be set explicitly for other than magnetic field.");
      return nullptr;
    }
  }

  // electromagnetic fields
  G4ElectroMagneticField* elMagField = nullptr;
  if (equation >= kEqElectroMagnetic && equation <= kEqEMfieldWithEDM)
  {
    elMagField = dynamic_cast<G4ElectroMagneticField*>(fG4Field);
    if (elMagField == nullptr)
    {
      G4Exception(
        "G4FieldSetup::CreateEquation:", "GeomFieldParameters0001",
        FatalErrorInArgument, "Incompatible field and equation.\n"
        "The field type must be set explicitly for other than magnetic field.");
      return nullptr;
    }
  }

  // electromagnetic fields
  switch (equation)
  {
    case kEqMagnetic:
      return new G4Mag_UsualEqRhs(magField);
      break;

    case kEqMagneticWithSpin:
      return new G4Mag_SpinEqRhs(magField);
      break;

    case kEqElectroMagnetic:
      return new G4EqMagElectricField(elMagField);
      break;

    case kEqEMfieldWithSpin:
      return new G4EqEMFieldWithSpin(elMagField);
      break;

    case kEqEMfieldWithEDM:
      return new G4EqEMFieldWithEDM(elMagField);
      break;

    // other fields
    case kEqGravity:
    case kEqMonopole:
    case kEqReplate:
      G4Exception(
        "G4FieldSetup::CreateEquation:", "GeomFieldParameters0001",
        FatalErrorInArgument, "Limitation: Equation not supported in G4FieldBuilder.\n"
        "Only magnetic and electromagnetic field can be constructed with G4FieldBuilder.");
      return nullptr;
      break;
    case kUserEquation:
      // nothing to be done
      return nullptr;
      break;
  }

  G4Exception(
    "G4FieldSetup::CreateEquation:", "GeomFieldParameters0001",
    FatalErrorInArgument, "Unknown equation type.");
  return nullptr;
}

//_____________________________________________________________________________
G4MagIntegratorStepper* G4FieldSetup::CreateStepper(
  G4EquationOfMotion* equation, G4StepperType stepper)
{
  // Set the integrator of particle's equation of motion

  // Check steppers which require equation of motion of G4Mag_EqRhs type
  auto eqRhs = dynamic_cast<G4Mag_EqRhs*>(equation);
  if ((eqRhs == nullptr) && (stepper > kTsitourasRK45))
  {
    G4Exception(
      "G4FieldSetup::CreateStepper:", "GeomFieldParameters0001",
      FatalErrorInArgument, 
      "The stepper type requires equation of motion of G4Mag_EqRhs type.");
    return nullptr;
  }

  switch (stepper)
  {
    case kBogackiShampine23:
      return new G4BogackiShampine23(equation);
      break;

    case kBogackiShampine45:
      return new G4BogackiShampine45(equation);
      break;

    case kCashKarpRKF45:
      return new G4CashKarpRKF45(equation);
      break;

    case kClassicalRK4:
      return new G4ClassicalRK4(equation);
      break;

    case kDoLoMcPriRK34:
      return new G4DoLoMcPriRK34(equation);
      break;

    case kDormandPrince745:
      return new G4DormandPrince745(equation);
      break;

    case kDormandPrinceRK56:
      return new G4DormandPrinceRK56(equation);
      break;

    case kDormandPrinceRK78:
      return new G4DormandPrinceRK78(equation);
      break;

    case kExplicitEuler:
      return new G4ExplicitEuler(equation);
      break;

    case kImplicitEuler:
      return new G4ImplicitEuler(equation);
      break;

    case kSimpleHeum:
      return new G4SimpleHeum(equation);
      break;

    case kSimpleRunge:
      return new G4SimpleRunge(equation);
      break;

    case kTsitourasRK45:
      return new G4TsitourasRK45(equation);
      break;

    case kConstRK4:
      return new G4ConstRK4(eqRhs);
      break;

    case kExactHelixStepper:
      return new G4ExactHelixStepper(eqRhs);
      break;

    case kHelixExplicitEuler:
      return new G4HelixExplicitEuler(eqRhs);
      break;

    case kHelixHeum:
      return new G4HelixHeum(eqRhs);
      break;

    case kHelixImplicitEuler:
      return new G4HelixImplicitEuler(eqRhs);
      break;

    case kHelixMixedStepper:
      return new G4HelixMixedStepper(eqRhs);
      break;

    case kHelixSimpleRunge:
      return new G4HelixSimpleRunge(eqRhs);
      break;

    case kNystromRK4:
      return new G4NystromRK4(eqRhs);
      break;

    case kRKG3Stepper:
      return new G4RKG3_Stepper(eqRhs);
      break;
    case kUserStepper:
      // nothing to be done
      return nullptr;
      break;

    // templated steppers
    case kTCashKarpRKF45:
    case kTDormandPrince45:
    case kTMagErrorStepper:
    case kQSStepper:
      G4Exception(
        "G4FieldSetup::CreateStepper:", "GeomFieldParameters0001",
        FatalErrorInArgument, "Limitation: Templated steppers not supported in G4FieldBuilder");
      return nullptr;
      break;

    default:
      G4Exception(
        "G4FieldSetup::CreateStepper:", "GeomFieldParameters0001",
        FatalErrorInArgument, "Incorrect stepper type.");
      return nullptr;
  }
}

//_____________________________________________________________________________
G4VIntegrationDriver* G4FieldSetup::CreateFSALStepperAndDriver(
  G4EquationOfMotion* equation, G4StepperType stepperType, G4double minStep)
{
  // Set the FSAL integrator of particle's equation of motion

  switch (stepperType)
  {
    case kRK547FEq1:
      return new G4FSALIntegrationDriver<G4RK547FEq1>(
        minStep, new G4RK547FEq1(equation));

    case kRK547FEq2:
      return new G4FSALIntegrationDriver<G4RK547FEq2>(
        minStep, new G4RK547FEq2(equation));

    case kRK547FEq3:
      return new G4FSALIntegrationDriver<G4RK547FEq3>(
        minStep, new G4RK547FEq3(equation));

    default:
      G4Exception(
        "G4FieldSetup::CreateFSALStepperAndDriver", "GeomFieldParameters0001",
        FatalErrorInArgument, "Incorrect stepper type.");
      return nullptr;
  }
}

//_____________________________________________________________________________
void G4FieldSetup::CreateCachedField()
{
  // Create cached field (if ConstDistance is set)
  fG4Field = CreateCachedField(fParameters, fG4Field);
}

//_____________________________________________________________________________
void G4FieldSetup::CreateStepper()
{
  // Create equation of motion (or get the user one if defined)
  if (fParameters.GetEquationType() == kUserEquation)
  {
    fEquation = fParameters.GetUserEquationOfMotion();
  }
  else
  {
    delete fEquation;
    fEquation = nullptr;
    fEquation = CreateEquation(fParameters.GetEquationType());
  }
  fEquation->SetFieldObj(fG4Field);

  // Create stepper  (or get the user one if defined)
  if (fParameters.GetStepperType() == kUserStepper)
  {
    // User stepper
    fStepper = fParameters.GetUserStepper();
  }
  else if (fParameters.GetStepperType() >= kRK547FEq1)
  {
    // FSAL stepper
    delete fDriver;
    delete fStepper;
    fDriver = nullptr;
    fStepper = nullptr;
    fDriver = CreateFSALStepperAndDriver(
      fEquation, fParameters.GetStepperType(), fParameters.GetMinimumStep());
    if (fDriver != nullptr)
    {
      fStepper = fDriver->GetStepper();
    }
  }
  else
  {
    // Normal stepper
    delete fStepper;
    fStepper = nullptr;
    fStepper = CreateStepper(fEquation, fParameters.GetStepperType());
  }
}

//_____________________________________________________________________________
void G4FieldSetup::CreateChordFinder()
{
  // Chord finder
  if (fParameters.GetFieldType() == kMagnetic)
  {
    if (fDriver != nullptr)
    {
      fChordFinder = new G4ChordFinder(fDriver);
    }
    else
    {
      // Chord finder
      fChordFinder = new G4ChordFinder(static_cast<G4MagneticField*>(fG4Field),
        fParameters.GetMinimumStep(), fStepper);
    }
    fChordFinder->SetDeltaChord(fParameters.GetDeltaChord());
  }
  else if (fParameters.GetFieldType() == kElectroMagnetic)
  {
    auto  intDriver = new G4MagInt_Driver(
      fParameters.GetMinimumStep(), fStepper, fStepper->GetNumberOfVariables());
    if (intDriver != nullptr)
    {
      // Chord finder
      fChordFinder = new G4ChordFinder(intDriver);
    }
  }
}

//_____________________________________________________________________________
void G4FieldSetup::UpdateFieldManager()
{
  fFieldManager->SetChordFinder(fChordFinder);
  fFieldManager->SetDetectorField(fG4Field);

  fFieldManager->SetMinimumEpsilonStep(fParameters.GetMinimumEpsilonStep());
  fFieldManager->SetMaximumEpsilonStep(fParameters.GetMaximumEpsilonStep());
  fFieldManager->SetDeltaOneStep(fParameters.GetDeltaOneStep());
  fFieldManager->SetDeltaIntersection(fParameters.GetDeltaIntersection());
}

//
// public methods
//

//_____________________________________________________________________________
void G4FieldSetup::Clear()
{
  // First clean up previous state.
  delete fChordFinder;
  fChordFinder = nullptr;

  if (fG4Field == nullptr)
  {
    delete fEquation;
    delete fDriver;
    delete fStepper;
    delete fChordFinder;
    fEquation = nullptr;
    fDriver = nullptr;
    fStepper = nullptr;
    fChordFinder = nullptr;
    fFieldManager->SetChordFinder(fChordFinder);
    fFieldManager->SetDetectorField(fG4Field);
  }
}

//_____________________________________________________________________________
void G4FieldSetup::Update()
{
  // Update field with new field parameters

  Clear();
  if (fG4Field == nullptr)
  {
    // No further update needed
    return;
  }

  CreateCachedField();
  CreateStepper();
  CreateChordFinder();
  UpdateFieldManager();
}

//_____________________________________________________________________________
void G4FieldSetup::PrintInfo(G4int verboseLevel, const G4String& about)
{
  if (verboseLevel == 0) { return; }

  auto fieldType = G4FieldParameters::FieldTypeName(fParameters.GetFieldType());
  auto isCachedMagneticField = (fParameters.GetConstDistance() > 0.);
  if (fLogicalVolume == nullptr)
  {
    fieldType = "Global";
  }
  else
  {
    fieldType = "Local (in ";
    fieldType.append(fLogicalVolume->GetName());
    fieldType.append(")");
  }
  if (isCachedMagneticField)
  {
    fieldType.append(" cached");
  }

  G4cout << fieldType << " field " << about << " with stepper ";
  G4cout << G4FieldParameters::StepperTypeName(fParameters.GetStepperType())
         << G4endl;

  if (verboseLevel > 1)
  {
    fParameters.PrintParameters();
  }
}
