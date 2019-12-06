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
//
/// \file field/field06/src/F06DetectorConstruction.cc
/// \brief Implementation of the F06DetectorConstruction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F06DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UserLimits.hh"
#include "G4SystemOfUnits.hh"

#include "G4UniformGravityField.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "G4RepleteEofM.hh"
//#include "G4EqGravityField.hh"

#include "G4ClassicalRK4.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4IntegrationDriver.hh"
#include "G4ChordFinder.hh"
#include "G4PropagatorInField.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F06DetectorConstruction::F06DetectorConstruction()
 : fVacuum(0), fSolidWorld(0), fLogicWorld(0), fPhysiWorld(0)
{
  // materials
  DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F06DetectorConstruction::~F06DetectorConstruction()
{
  if (fField) delete fField;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F06DetectorConstruction::DefineMaterials()
{
  G4NistManager* nistMan = G4NistManager::Instance();

  fVacuum = nistMan->FindOrBuildMaterial("G4_Galactic");

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* F06DetectorConstruction::Construct()
{
  //
  // World
  //

  G4double expHall_x = 1.0*m;
  G4double expHall_y = 1.0*m;
  G4double expHall_z = 1.0*m;

  fSolidWorld = new G4Box("World",                 //its name
                   expHall_x,expHall_y,expHall_z); //its size

  fLogicWorld = new G4LogicalVolume(fSolidWorld,   //its solid
                                    fVacuum,       //its material
                                    "World");      //its name

  fPhysiWorld = new G4PVPlacement(0,               //no rotation
                                  G4ThreeVector(), //at (0,0,0)
                                  fLogicWorld,     //its logical volume
                                  "World",         //its name
                                  0,               //its mother  volume
                                  false,           //no boolean operation
                                  0);              //copy number

  G4double maxStep = 1.0*mm;
  G4double maxTime = 41.*s;

  G4UserLimits* stepLimit = new G4UserLimits(maxStep,DBL_MAX,maxTime);

  fLogicWorld->SetUserLimits(stepLimit);
 
  //
  // Visualization attributes
  //
  // fLogicWorld->SetVisAttributes (G4VisAttributes::GetInvisible());

  //
  //always return the physical World
  //
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4UniformGravityField* F06DetectorConstruction::fField = 0;

void F06DetectorConstruction::ConstructSDandField()
{
  using StepperType = G4ClassicalRK4;
   
  if (!fField) {

     fField = new G4UniformGravityField();

     G4RepleteEofM* equation = new G4RepleteEofM(fField);
//     G4EqGravityField* equation = new G4EqGravityField(fField);

     G4TransportationManager* transportMgr =
                           G4TransportationManager::GetTransportationManager();

     G4FieldManager* fieldManager= transportMgr->GetFieldManager();
     fieldManager->SetDetectorField(fField);

     const int nVar= 8;  // 12 for RepleteEofM
     StepperType* stepper = new StepperType(equation,nVar);

     G4double minStep           = 0.01*mm;
     G4ChordFinder* chordFinder = nullptr;     
     if( stepper )
     {
        auto intgrDriver = new G4IntegrationDriver<StepperType>(minStep,
                                          stepper,
                                          stepper->GetNumberOfVariables());
        if( intgrDriver ){
           chordFinder = new G4ChordFinder(intgrDriver);
        }
     }

     // OLD -- and wrong
     // new G4ChordFinder((G4MagneticField*)fField,minStep,stepper);
     
     // Set accuracy parameters
     G4double deltaChord        = 3.0*mm;
     chordFinder->SetDeltaChord( deltaChord );

     G4double deltaIntersection = 0.1*mm;
     fieldManager->SetDeltaIntersection(deltaIntersection);

     //  Control accuracy of integration
     //
     G4double deltaOneStep = 0.01*mm;
     fieldManager->SetAccuraciesWithDeltaOneStep(deltaOneStep);
     //  
     G4double epsMax       = 1.0e-4;  // Pure number -- maximum relative integration error
     G4double epsMin       = 2.5e-7;  // 
     fieldManager->SetMinimumEpsilonStep(epsMin);
     fieldManager->SetMaximumEpsilonStep(epsMax);
     // The acceptable relative accuracy is calculated  as  deltaOneStep / stepsize
     //    but bounded to the interval between these values!
     
     fieldManager->SetChordFinder(chordFinder);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
