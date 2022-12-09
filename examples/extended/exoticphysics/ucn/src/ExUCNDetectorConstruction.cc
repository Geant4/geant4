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
/// \file exoticphysics/ucn/src/ExUCNDetectorConstruction.cc
/// \brief Implementation of the ExUCNDetectorConstruction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExUCNDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4UCNMaterialPropertiesTable.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
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
#include "G4PhysicalConstants.hh"

#include "G4UniformGravityField.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "G4RepleteEofM.hh"
//#include "G4EqGravityField.hh"

#include "G4ClassicalRK4.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"
#include "G4PropagatorInField.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExUCNDetectorConstruction::ExUCNDetectorConstruction()
 : fVacuum(0), fGuideMaterial(0)
{
  // materials
  DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExUCNDetectorConstruction::~ExUCNDetectorConstruction()
{
  if (fField) delete fField;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExUCNDetectorConstruction::DefineMaterials()
{
  G4NistManager* nistMan = G4NistManager::Instance();

  fVacuum = nistMan->FindOrBuildMaterial("G4_Galactic");
  fGuideMaterial = nistMan->FindOrBuildMaterial("G4_Ni");

  // --- Ni diffuse 10%

  G4UCNMaterialPropertiesTable* MPT = new G4UCNMaterialPropertiesTable();

  //  MPT->AddConstProperty("REFLECTIVITY",1.); 
  //  Commented out above line as REFLECTIVITY=1 by default in
  //  G4OpBoundaryProcess.  Also use AddProperty to set REFLECTIVITY if needed
  MPT->AddConstProperty("DIFFUSION",0.1);
  MPT->AddConstProperty("FERMIPOT",252.0); // Gollub, Table 2.1 in neV
  MPT->AddConstProperty("SPINFLIP",0.);
  MPT->AddConstProperty("LOSS", 12.5e-5); //  Gollub, Table 2.1
  MPT->AddConstProperty("LOSSCS",0.);
  MPT->AddConstProperty("ABSCS",4.49); // 1/v loss cross-section  at room temp.
  MPT->AddConstProperty("SCATCS",18.5); // (incoherent) "elastic" scattering cs

  G4double neV = 1.e-9*eV;

  MPT->SetMicroRoughnessParameters(30*nm, 1*nm,
                                   180, 1000,
                                   0*degree, 90*degree,
                                   1*neV, 1000*neV,
                                   15, 15,
                                   0.01*degree);

  fGuideMaterial->SetMaterialPropertiesTable(MPT);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ExUCNDetectorConstruction::Construct()
{
  //
  // World
  //

  G4double worldSizeX =   1.*m;
  G4double worldSizeY =   1.*m;
  G4double worldSizeZ = 100.*m;

  G4Box* solidWorld = new G4Box("World",
                                worldSizeX/2.,worldSizeY/2.,worldSizeZ/2.);

  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,
                                                    fVacuum,
                                                    "World");

  G4VPhysicalVolume* physiWorld = new G4PVPlacement(0,
                                                    G4ThreeVector(),
                                                    "World",
                                                    logicWorld,
                                                    0,
                                                    false,
                                                    0);

// --------------------------------- Guide -------------------------------------

  G4double GuideR = 35.*mm;
  G4double GuideW =  2.*mm;
  G4double GuideL =  6.*m;

  G4Tubs* solidGuide = new G4Tubs("SolidGuide",
                                  GuideR,GuideR+GuideW,GuideL/2.,0.,twopi);

  G4LogicalVolume* logicGuide = new G4LogicalVolume(solidGuide,
                                                    fGuideMaterial,
                                                    "Guide");

  new G4PVPlacement(0,G4ThreeVector(),"Guide",logicGuide,physiWorld,false,0);

// ------------------------------ End Plate  -----------------------------------

  G4Tubs* solidEndPlate = new G4Tubs("EndPlate",0.,GuideR,GuideW/2.,0.,twopi);

  G4LogicalVolume* logicEndPlate = new G4LogicalVolume(solidEndPlate,
                                                       fVacuum,
                                                       "EndPlate");

  G4ThreeVector endPlatePos = G4ThreeVector(0.,0.,GuideL/2.+GuideW/2.);

  new G4PVPlacement(0,endPlatePos,"EndPlate",logicEndPlate,physiWorld,false,0);

  G4double maxStep = 1.0*mm;
  G4double maxTime = 100.*s;

  G4UserLimits* stepLimit = new G4UserLimits(maxStep,DBL_MAX,maxTime);

  logicWorld->SetUserLimits(stepLimit);
 
  //
  // Visualization attributes
  //

  G4VisAttributes* guideColor = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  guideColor->SetVisibility(true);
  guideColor->SetForceWireframe(true);

  G4VisAttributes* endPlateColor = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  endPlateColor->SetVisibility(true);
  endPlateColor->SetForceSolid(true);

  logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());
  logicGuide->SetVisAttributes(guideColor);
  logicEndPlate->SetVisAttributes(endPlateColor);

  //
  //always return the physical World
  //
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4UniformGravityField* ExUCNDetectorConstruction::fField = 0;

void ExUCNDetectorConstruction::ConstructSDandField()
{
  if (!fField) {

     fField = new G4UniformGravityField();

     G4RepleteEofM* equation = new G4RepleteEofM(fField);
//     G4RepleteEofM* equation = new G4RepleteEofM(fField,12);
//     G4EqGravityField* equation = new G4EqGravityField(fField);

     G4FieldManager* fieldManager
      = G4TransportationManager::GetTransportationManager()->GetFieldManager();
     fieldManager->SetDetectorField(fField);

     G4MagIntegratorStepper* stepper = new G4ClassicalRK4(equation,8);
//     G4MagIntegratorStepper* stepper = new G4ClassicalRK4(equation,12);

     G4double minStep           = 0.01*mm;

     G4ChordFinder* chordFinder = 
                   new G4ChordFinder((G4MagneticField*)fField,minStep,stepper);

     // Set accuracy parameters
     G4double deltaChord        = 3.0*mm;
     chordFinder->SetDeltaChord( deltaChord );

     G4double deltaOneStep      = 0.01*mm;
     fieldManager->SetAccuraciesWithDeltaOneStep(deltaOneStep);

     G4double deltaIntersection = 0.1*mm;
     fieldManager->SetDeltaIntersection(deltaIntersection);

     G4TransportationManager* transportManager =
                           G4TransportationManager::GetTransportationManager();

     G4PropagatorInField* fieldPropagator =
                                      transportManager->GetPropagatorInField();

     // Dimensionless limits for relative accuracy of integration
     G4double epsMin            = 2.5e-7;
     G4double epsMax            = 0.001; // Will soon be maximum without warning.
     // The relative accuracy used for a step of length 'l'
     //                    a.)  epsMin              if deltaOneStep / l < epsMin
     //    epsilon_step =  b.)  epsMax              if deltaOneStep / l > epsMax
     //                    c.)  deltaOneStep / l    otherwise

     fieldPropagator->SetMinimumEpsilonStep(epsMin);
     fieldPropagator->SetMaximumEpsilonStep(epsMax);

     fieldManager->SetChordFinder(chordFinder);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
