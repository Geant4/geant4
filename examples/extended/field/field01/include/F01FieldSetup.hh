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
/// \file field/field01/include/F01FieldSetup.hh
/// \brief Definition of the F01FieldSetup class
//
//
//
//
//  A class for control of the Magnetic Field of the detector.
//  The field is assumed to be uniform.
//
// Should this be a:
//    i) a messenger class
//   ii) user class that creates the field           ( Current choice )
//  iii) a field manager that creates/updates field
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef F01FieldSetup_h
#define F01FieldSetup_h 1

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"

class G4FieldManager;
class G4ChordFinder;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
class G4VIntegrationDriver;
class F01FieldMessenger;

class F01FieldSetup
{
public:
  F01FieldSetup(G4ThreeVector,                //  The value of the field
                G4int  stepperNum = -1000,    //  -ive = fsal, +ive = old.
                                              // default (-1000) uses next flag
                G4bool useFSALstepper= false );

  F01FieldSetup();               //  A zero field
  F01FieldSetup( F01FieldSetup & ) = delete;

  virtual ~F01FieldSetup();

  void SetStepperType( G4int i )
     { fStepperType = i; CreateStepperAndChordFinder(); }

  void SetStepper();

  void SetMinStep(G4double s) { fMinStep = s; }

  void InitialiseAll();    //  Set parameters and call method below

   // Original method
  void CreateStepperAndChordFinder();

   // Create FSAL stepper and driver
  void CreateFSALStepperAndChordFinder();

   // Create Boris driver ( 2nd order Symplectic integration )
  void CreateAndSetupBorisDriver();

  void   SetUseFSALstepper(G4bool val= true) { fUseFSALstepper = val; }
  G4bool GetUseFSALstepper()                 { return fUseFSALstepper; }

  // Parameters for integration accuracy : get / control
  // -------------------------------------------------------
  G4double GetEpsilonMin(){ return fDesiredEpsilonMin; }
  G4double GetEpsilonMax(){ return fDesiredEpsilonMax; }

  void   SetEpsilonMin(G4double val){ fDesiredEpsilonMin= val; }
  void   SetEpsilonMax(G4double val){ fDesiredEpsilonMax= val; }

  void   SetDeltaOneStep(G4double val){ fDeltaOneStep= val; }

  // Special for this setup only
  // ----------------------------
  void SetFieldValue(G4ThreeVector fieldVector);
  void SetFieldZValue(G4double      fieldValue);
  G4ThreeVector GetConstantFieldValue();

protected:
   // Implementation methods
  G4VIntegrationDriver* CreateFSALStepperAndDriver();

   // Find the global Field Manager
  G4FieldManager*          GetGlobalFieldManager();

protected:

  G4FieldManager*          fFieldManager = nullptr;
  G4ChordFinder*           fChordFinder = nullptr;
  G4Mag_UsualEqRhs*        fEquation = nullptr;
  G4MagneticField*         fMagneticField = nullptr;

  G4MagIntegratorStepper*  fStepper = nullptr;
  G4bool                   fUseFSALstepper = false;
  G4VIntegrationDriver*    fDriver =  nullptr;  // If non-null, its new type (FSAL)
  G4int                    fStepperType = -1;
  F01FieldMessenger*       fFieldMessenger = nullptr;

  // Parameters / Invariant during tracking loop
  G4double             fMinStep     = -1.0;
  G4double             fDeltaOneStep= -1.0;
  G4double             fDesiredEpsilonMin = 1.0e-05; // tight: 1.0e-8  std: 1.0e-5 to 1.e-6  loose: 1.0e-4
  G4double             fDesiredEpsilonMax = 0.005;   // tight: 1.0e-5  std: 1.0e-4 to 5.e-3  loose: 1.0e-3+
};

#endif
