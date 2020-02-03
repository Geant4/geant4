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
/// \file electromagnetic/TestEm10/src/DetectorALICE06.cc
/// \brief Implementation of the DetectorALICE06 class
//
//
//
//

#include "DetectorALICE06.hh"
#include "SensitiveDetector.hh"
#include "Materials.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"

#include "G4Region.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorALICE06::DetectorALICE06()
  : fRadiatorDescription(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorALICE06::~DetectorALICE06()
{
  // delete fRadiatorDescription;
        // the description is deleted in detector construction
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorALICE06::Construct()
{
  // Geometry parameters
  //

  G4cout << "DetectorALICE06 setup" << G4endl;

  G4double worldSizeZ = 600.*cm;
  G4double worldSizeR = 22.*cm;

  // Radiator and detector parameters

  G4double radThickness = 0.020*mm;
  G4double gasGap       = 0.500*mm;  
  G4double foilGasRatio = radThickness/(radThickness+gasGap);
  G4int foilNumber      = 120;

  G4double absorberThickness = 37*mm;
  G4double absorberRadius    = 100.*mm;

  G4double electrodeThick = 100.0*micrometer;
  G4double pipeLength     = 160.0*cm;
  G4double mylarThick     = 20.0*micrometer;
  G4double detGap         = 0.01*mm;

  G4double startZ = 100.0*mm;

  // Materials
  //

  // Change to create materials using NIST
  G4Material* air   = Materials::GetInstance()->GetMaterial("Air");
  G4Material* ch2   = Materials::GetInstance()->GetMaterial("CH2");
  G4Material* xe15CO2 = Materials::GetInstance()->GetMaterial("Xe15CO2");

  G4double foilDensity = ch2->GetDensity();
  G4double gasDensity  = air->GetDensity();  
  G4double totDensity  = foilDensity*foilGasRatio 
                       + gasDensity*(1.0-foilGasRatio);

  G4double fractionFoil =  foilDensity*foilGasRatio/totDensity;
  G4double fractionGas  =  1.0 - fractionFoil;
  G4Material* radiatorMat = new G4Material("radiatorMat", totDensity, 2);
  radiatorMat->AddMaterial(ch2, fractionFoil);
  radiatorMat->AddMaterial(air, fractionGas);

  // Radiator description
  fRadiatorDescription = new RadiatorDescription;
  fRadiatorDescription->fFoilMaterial  = ch2;   // CH2; // Kapton; // Mylar ; // Li ; // CH2 ;
  fRadiatorDescription->fGasMaterial   = air;   // CO2; // He; //
  fRadiatorDescription->fFoilThickness = radThickness;
  fRadiatorDescription->fGasThickness  = gasGap;
  fRadiatorDescription->fFoilNumber = foilNumber;
  
  G4Material* worldMaterial    = air; // CO2;
  G4Material* absorberMaterial = xe15CO2;

  // Volumes
  //
 
  G4VSolid* solidWorld 
    = new G4Box("World", worldSizeR, worldSizeR, worldSizeZ/2.);
 
  G4LogicalVolume* logicWorld 
    = new G4LogicalVolume(solidWorld,  worldMaterial,  "World");

  G4VPhysicalVolume* physicsWorld 
    = new G4PVPlacement(0, G4ThreeVector(), "World", logicWorld, 0,  false, 0);

  // TR radiator envelope

  G4double radThick = foilNumber*(radThickness + gasGap) - gasGap + detGap;
  G4double radZ = startZ + 0.5*radThick;

  G4VSolid* solidRadiator 
    = new G4Box("Radiator", 1.1*absorberRadius, 1.1*absorberRadius, 0.5*radThick);

  G4LogicalVolume* logicRadiator 
    = new G4LogicalVolume(solidRadiator, radiatorMat, "Radiator");
 
  new G4PVPlacement(0, G4ThreeVector(0, 0, radZ),
                    "Radiator", logicRadiator, physicsWorld, false, 0 );
  
  fRadiatorDescription->fLogicalVolume = logicRadiator;

  // Create region for radiator

  G4Region* radRegion = new G4Region("XTRradiator");
  radRegion->AddRootLogicalVolume(logicRadiator);

  // Drift Electrode on both sides of Radiator
  // (not placed)

  G4double zElectrode1 = radZ - radThick/2. - electrodeThick/2.;
  G4double zElectrode2 = radZ + radThick/2. + electrodeThick/2.;

  G4cout << "zElectrode1 = " << zElectrode1/mm << " mm" <<G4endl;
  G4cout << "zElectrode2 = " << zElectrode2/mm << " mm" <<G4endl;
  G4cout << "fElectrodeThick = " << electrodeThick/mm << " mm" << G4endl <<G4endl;

  // Helium Pipe
  // (not placed)

  //Distance between pipe and radiator / absorber
  G4double pipeDist      = 1.*cm;
  G4double zPipe = zElectrode2 + electrodeThick/2. + pipeDist/2. + pipeLength/2.;

  G4cout << "zPipe = " << zPipe/mm << " mm" << G4endl;
  G4cout << "pipeLength = " << pipeLength/mm <<" mm" << G4endl << G4endl;

  // Mylar Foil on both sides of helium pipe
  // (not placed)

  G4double zMylar1 = zPipe - pipeLength/2. - mylarThick/2. - 0.001*mm;
  G4double zMylar2 = zPipe + pipeLength/2. + mylarThick/2. + 0.001*mm;

  G4cout << "zMylar1 = " << zMylar1/mm << " mm" << G4endl;
  G4cout << "zMylar2 = " << zMylar2/mm << " mm" << G4endl;
  G4cout << "fMylarThick = " << mylarThick/mm << " mm" << G4endl << G4endl;

  // Mylar Foil on Chamber
  // (not placed)

  G4double zMylar = zElectrode2 + electrodeThick/2. + mylarThick/2. + 1.0*mm;
  zMylar += ( pipeLength + pipeDist );

  G4cout << "zMylar = " << zMylar/mm <<" mm" <<G4endl;
  G4cout << "mylarThick = " << mylarThick/mm << " mm" << G4endl << G4endl;

  // Absorber

  G4double absorberZ = zMylar + mylarThick + absorberThickness/2.;

  G4VSolid* solidAbsorber 
    = new G4Box("Absorber", absorberRadius, 10.*mm, absorberThickness/2.);

  G4LogicalVolume* logicAbsorber 
    = new G4LogicalVolume(solidAbsorber, absorberMaterial, "Absorber");

  new G4PVPlacement(0, G4ThreeVector(0., 0., absorberZ),
                    "Absorber", logicAbsorber, physicsWorld, false, 0);

  G4Region* regGasDet = new G4Region("XTRdEdxDetector");
  regGasDet->AddRootLogicalVolume(logicAbsorber);

  // Sensitive Detectors: Absorber

  SensitiveDetector* sd = new SensitiveDetector("AbsorberSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(sd );
  logicAbsorber->SetSensitiveDetector(sd);

  // Print geometry parameters

  G4cout << "\n The  WORLD   is made of "
         << worldSizeZ/mm << "mm of " << worldMaterial->GetName();
  G4cout << ", the transverse size (R) of the world is " 
         << worldSizeR/mm << " mm. " << G4endl;
  G4cout << " The ABSORBER is made of "
         << absorberThickness/mm << "mm of " << absorberMaterial->GetName();
  G4cout << ", the transverse size (R) is " 
         << absorberRadius/mm << " mm. " << G4endl;
  G4cout << " Z position of the (middle of the) absorber " 
         << absorberZ/mm << "  mm." << G4endl;

  G4cout << "radZ = " << radZ/mm << " mm" << G4endl;
  G4cout << "startZ = " << startZ/mm<< " mm" << G4endl;

  G4cout << "fRadThick = " << radThick/mm << " mm"<<G4endl;
  G4cout << "fFoilNumber = " << foilNumber << G4endl;
  G4cout << "fRadiatorMat = " << radiatorMat->GetName() << G4endl;
  G4cout << "WorldMaterial = " << worldMaterial->GetName() << G4endl;
  G4cout << G4endl;

  return physicsWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
