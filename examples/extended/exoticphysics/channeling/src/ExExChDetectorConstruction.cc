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
/// \file channeling/src/ExExChDetectorConstruction.cc
/// \brief Implementation of the ExExChDetectorConstruction class
//
//

#include "ExExChDetectorConstruction.hh"
#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SDManager.hh"
#include "ExExChSensitiveDetector.hh"

#include "XLatticeManager3.hh"

#include "XLogicalAtomicLattice.hh"
#include "XLogicalAtomicLatticeDiamond.hh"
#include "XLogicalBase.hh"
#include "XUnitCell.hh"

#include "XCrystalPlanarMolierePotential.hh"
#include "XCrystalPlanarMoliereElectricField.hh"
#include "XCrystalPlanarNucleiDensity.hh"
#include "XCrystalPlanarMoliereElectronDensity.hh"
#include "XCrystalPlanarMoliereTempPotential.hh"
#include "XCrystalPlanarMoliereElectricField.hh"
#include "XCrystalPlanarNucleiDensity.hh"
#include "XCrystalPlanarMoliereElectronDensity.hh"
#include "XCrystalCharacteristicArray.hh"
#include "XCrystalIntegratedDensityPlanar.hh"
#include "XCrystalIntegratedDensityHub.hh"
#include "XCrystalIntegratedDensityHub.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExExChDetectorConstruction::ExExChDetectorConstruction():
fWorldLogic(0),fXtalLogic(0){
    
    //SiSD standard parameters
    fSiSD = true;
    fSSDSize = G4ThreeVector(3.8 * CLHEP::centimeter,
                             3.8 * CLHEP::centimeter,
                             640. * CLHEP::micrometer);
    fSSDXtalDistance[0] = - (9998.) * CLHEP::millimeter;
    fSSDXtalDistance[1] = - (320.) * CLHEP::millimeter;
    fSSDXtalDistance[2] = + (10756.) * CLHEP::millimeter;

    //SiSD Box standard parameters
    fSSDBoxSize = G4ThreeVector(25. * CLHEP::centimeter,
                                25. * CLHEP::centimeter,
                                10. * CLHEP::centimeter);
    fSSDBoxThickness = 4. * CLHEP::millimeter;

    //Beampipe standard parameters
    fBeamPipe = false;
    fBeamPipeThickness = (0.3) * CLHEP::centimeter;
    fBeamPipeRadius = (15.6) * CLHEP::centimeter;
    fXtal = true;
    fXtalAngle = G4ThreeVector(0.,0.,0.);
    fXtalSize = G4ThreeVector(1. * CLHEP::millimeter,
                              70. * CLHEP::millimeter,
                              1.94 * CLHEP::millimeter);
    fXtalCurvatureRadius = G4ThreeVector(38.416 * CLHEP::meter,
                                         0. * CLHEP::meter,
                                         0. * CLHEP::meter);
    
    fXtalCellSize = G4ThreeVector(5.431 * CLHEP::angstrom,
                                  5.431 * CLHEP::angstrom,
                                  5.431 * CLHEP::angstrom);
    
    fXtalCellAngle = G4ThreeVector(90.*CLHEP::deg,
                                   90.*CLHEP::deg,
                                   90.*CLHEP::deg);
    
    fXtalTVA = 0.075 * CLHEP::angstrom;
    fXtalMiller = G4ThreeVector(2,2,0);
    
    SetXtalMaterial("G4_Si");
    
    fMessenger = new ExExChDetectorConstructionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ExExChDetectorConstruction::~ExExChDetectorConstruction(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ExExChDetectorConstruction::DefineMaterials(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ExExChDetectorConstruction::Construct(){
    //** World **//
    fWorldSize = G4ThreeVector(1. * CLHEP::meter,
                               1. * CLHEP::meter,
                               30. * CLHEP::meter);
    fWorldMaterial = G4NistManager::
        Instance()->FindOrBuildMaterial("G4_Galactic");
    
    fWorldSolid = new G4Box("World",
                            fWorldSize.x()/2.,
                            fWorldSize.y()/2.,
                            fWorldSize.z()/2.);
    
    fWorldLogic = new G4LogicalVolume(fWorldSolid,
                                      fWorldMaterial,
                                      "World");
    
    fWorldPhysical = new G4PVPlacement(0,
                                       G4ThreeVector(),
                                       fWorldLogic,
                                       "World",
                                       0,
                                       false,
                                       0);
    
    //** SiSD **//

    if(fSiSD){
        for(unsigned int i1=0;i1<3;i1++){
                G4LogicalVolume* fSSDBoxLogic = ConstructSiSD(i1);

                G4ThreeVector vBoxPosition = 
                G4ThreeVector(+fSSDBoxSize.x()/4.,
                -fSSDBoxSize.y()/4.,
                fSSDXtalDistance[i1]);

            new G4PVPlacement(0,
                              vBoxPosition,
                              fSSDBoxLogic,"SiSD",
                              fWorldLogic,
                              false,
                              i1);
        }
    }
    //** BeamPipe **//
    if(fBeamPipe){
        G4double fBeamPipeFromSiSDDistance = 20. * CLHEP::centimeter;

        G4ThreeVector fBeamPipeA0Position =
                G4ThreeVector(0.,
                              0.,
                              fSSDXtalDistance[0] + std::fabs(fSSDXtalDistance[1] - 
                              fSSDXtalDistance[0])/2.);
                              
        G4double fBeamPipeA0Length =
                std::fabs(fSSDXtalDistance[1] - fSSDXtalDistance[0]) -
                2. * (fSSDSize.z()/2.) - 2. * fBeamPipeFromSiSDDistance;

            G4LogicalVolume* fBeamPipeA0Logic =
                    ConstructBeamPipe(fBeamPipeA0Length);

        new G4PVPlacement(0,
                          fBeamPipeA0Position,
                          fBeamPipeA0Logic,
                         "BeamPipeA0",
                         fWorldLogic,
                         false,
                         0);

        G4ThreeVector fBeamPipeA1Position =
                G4ThreeVector(0.,0.,+ std::fabs(fSSDXtalDistance[2]) /2.);
        G4double fBeamPipeA1Length =
                std::fabs(fSSDXtalDistance[2]) - 2. * (fSSDSize.z()/2.) -
                2. * fBeamPipeFromSiSDDistance;
                
            G4LogicalVolume* fBeamPipeA1Logic =
                    ConstructBeamPipe(fBeamPipeA1Length);

        new G4PVPlacement(0,
                          fBeamPipeA1Position,
                          fBeamPipeA1Logic,
                          "BeamPipeA1",
                          fWorldLogic,
                          false,
                          1);
    }
#ifndef G4MULTITHREADED
    G4String SDname;
    G4VSensitiveDetector* telescope =
        new ExExChSensitiveDetector(SDname="/telescope");
    G4SDManager::GetSDMpointer()->AddNewDetector(telescope);
    for(unsigned int i1=0;i1<3;i1++){
            fSSDLogic[i1]->SetSensitiveDetector(telescope);
    }
#endif
    //** Crystal **//
    if(fXtal){
            ConstructXtalTarget();
    }
    
    return fWorldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef G4MULTITHREADED
void ExExChDetectorConstruction::ConstructSDandField(){
    G4String SDname;
    G4VSensitiveDetector* telescope =
        new ExExChSensitiveDetector(SDname="/telescope");
    G4SDManager::GetSDMpointer()->AddNewDetector(telescope);
    for(unsigned int i1=0;i1<3;i1++){
            fSSDLogic[i1]->SetSensitiveDetector(telescope);
    }
}
#else
void ExExChDetectorConstruction::ConstructSDandField(){
}
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* ExExChDetectorConstruction::ConstructSiSD(G4int copyNo){
    G4Material* Si = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
    G4Material* Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    G4Material* Galactic =
            G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    G4double fFoilThickness = 0.024 * CLHEP::millimeter;

    //** SiSD Box **//
    G4Box* fSSDBoxEnvelopeSolid =
            new G4Box("SiSDBoxEnvelope",
                  fSSDBoxSize.x()/2.,
                  fSSDBoxSize.y()/2.,
                  fSSDBoxSize.z()/2. + fFoilThickness);
    
    G4LogicalVolume* fSSDBoxEnvelopeLogic =
        new G4LogicalVolume(fSSDBoxEnvelopeSolid,
                            Galactic,
                            "SiSDBoxEnvelope");

    //** SiSD Active Area **//
    G4Box* fSSDSolid = new G4Box("SiSD",
                          fSSDSize.x()/2.,
                          fSSDSize.y()/2.,
                          fSSDSize.z()/2.);
    fSSDLogic[copyNo] = new G4LogicalVolume(fSSDSolid,Si,"SiSD");

    //** SiSD Box **//
    G4Box* fSSDBoxSolidA = new G4Box("SiSDBoxA",
                                     fSSDBoxSize.x()/2.,
                                         fSSDBoxSize.y()/2.,
                                         fSSDBoxSize.z()/2.);

    G4Box* fSSDBoxSolidB = new G4Box("SiSDBoxB",
                    fSSDBoxSize.x()/2. - fSSDBoxThickness/2.,
                    fSSDBoxSize.y()/2. - fSSDBoxThickness/2.,
                    fSSDBoxSize.z()/2. - fSSDBoxThickness/2.);

    G4Box* fSSDBoxSolidC = new G4Box("SiSDBoxC",
                                         fSSDSize.x()/2. + fSSDBoxThickness/2.,
                                     fSSDSize.y()/2. + fSSDBoxThickness/2.,
                                     fSSDBoxSize.z()/2.);

    G4Box* fSSDBoxSolidD = new G4Box("SiSDBoxD",
                                     fSSDSize.x()/2.,
                                     fSSDSize.y()/2.,
                                     fSSDBoxSize.z());

    G4Box* fSSDBoxFoilSolid = new G4Box("SSDBoxFoil",
                                        fSSDSize.x()/2.,
                                        fSSDSize.y()/2.,
                                        fFoilThickness);
    G4SubtractionSolid* fSSDBoxSolid =
            new G4SubtractionSolid("SiSDBox",
                                   fSSDBoxSolidA,
                                   fSSDBoxSolidB);
                                           
    fSSDBoxSolid = new G4SubtractionSolid("SiSDBox",
                                          fSSDBoxSolid,
                                          fSSDBoxSolidC,
                                          0,
                                          G4ThreeVector(-fSSDBoxSize.x()/4.,
                                                        fSSDBoxSize.y()/4.,
                                                        0.));

    G4SubtractionSolid* fSSDBoxInternalSolid =
            new G4SubtractionSolid("SiSDBoxInternalSolid",
                               fSSDBoxSolidC,
                               fSSDBoxSolidD);

    G4LogicalVolume* fSSDBoxLogic = new G4LogicalVolume(fSSDBoxSolid,
                                                        Al,
                                                        "SiSDBox");
    G4LogicalVolume* fSSDBoxInternalLogic =
            new G4LogicalVolume(fSSDBoxInternalSolid,
                             Al,
                            "SiSDBox");
                            
    G4LogicalVolume* fSSDBoxFoilLogic = new G4LogicalVolume(fSSDBoxFoilSolid,
                                                            Al,
                                                            "SiSDBoxFoil");

    G4VisAttributes* fSSDBoxVisAttribute =
            new G4VisAttributes(G4Colour(0.7,0.7,0.7));          
    fSSDBoxVisAttribute->SetForceSolid(true);
    fSSDBoxLogic->SetVisAttributes(fSSDBoxVisAttribute);
    fSSDBoxInternalLogic->SetVisAttributes(fSSDBoxVisAttribute);

    G4VisAttributes* fSSDBoxFoilVisAttribute =
            new G4VisAttributes(G4Colour(0.8,0.8,0.8));
    fSSDBoxFoilVisAttribute->SetForceSolid(false);
    fSSDBoxFoilLogic->SetVisAttributes(fSSDBoxFoilVisAttribute);

    G4VisAttributes* fSiSDVisAttribute =
            new G4VisAttributes(G4Colour(1.0,0.65,0.0));
    fSiSDVisAttribute->SetForceSolid(true);
    fSSDLogic[copyNo]->SetVisAttributes(fSiSDVisAttribute);

    //** Add to Physical World **//
    new G4PVPlacement(0,
                      G4ThreeVector(-fSSDBoxSize.x()/4.,fSSDBoxSize.y()/4.,0.),
                      fSSDLogic[copyNo],"SiSD",
                      fSSDBoxEnvelopeLogic,
                      false,
                      copyNo);
     new G4PVPlacement(0,
                       G4ThreeVector(),
                       fSSDBoxLogic,"SiSDBox",
                       fSSDBoxEnvelopeLogic,
                       false,
                       copyNo);

     new G4PVPlacement(0,
                       G4ThreeVector(-fSSDBoxSize.x()/4.,fSSDBoxSize.y()/4.,0.),
                       fSSDBoxInternalLogic,"SiSDBox",
                       fSSDBoxEnvelopeLogic,
                       false,
                       copyNo);

     new G4PVPlacement(0,
                       G4ThreeVector(-fSSDBoxSize.x()/4.,
                                     fSSDBoxSize.y()/4.,
                                     (fSSDBoxSize.z()/2. - fFoilThickness/2.)),
                       fSSDBoxFoilLogic,"SiSDBoxFoil",
                       fSSDBoxEnvelopeLogic,
                       false,
                       copyNo);

     new G4PVPlacement(0,
                       G4ThreeVector(-fSSDBoxSize.x()/4.,
                                     fSSDBoxSize.y()/4.,
                                     -(fSSDBoxSize.z()/2. - fFoilThickness/2.)),
                       fSSDBoxFoilLogic,"SiSDBoxFoil",
                       fSSDBoxEnvelopeLogic,
                       false,
                       G4int(copyNo*2+1));

     return fSSDBoxEnvelopeLogic;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* ExExChDetectorConstruction::ConstructBeamPipe(G4double length){
    G4Material* Galactic =
            G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    G4Material* Mylar =
            G4NistManager::Instance()->FindOrBuildMaterial("G4_MYLAR");

    // Definition of vacuum
    G4double z = 7.;
    G4double a = 14.007*CLHEP::g/CLHEP::mole;
    G4double density = CLHEP::universe_mean_density;
    G4double pressure = 1.E-8 * 1.E-3 * CLHEP::bar;
    G4double temperature = 300.*CLHEP::kelvin;

    G4Material* Vacuum = new G4Material("Vacuum",
                                        z,
                                        a,
                                        density,
                                        kStateGas,
                                        temperature,
                                        pressure);

    // Definition of stainless steel (not in NIST) for pipes
    G4Element* elC = G4NistManager::Instance()->FindOrBuildElement("C");
    G4Element* elSi = G4NistManager::Instance()->FindOrBuildElement("Si");
    G4Element* elCr = G4NistManager::Instance()->FindOrBuildElement("Cr");
    G4Element* elMn = G4NistManager::Instance()->FindOrBuildElement("Mn");
    G4Element* elNi = G4NistManager::Instance()->FindOrBuildElement("Ni");
    G4Element* elFe = G4NistManager::Instance()->FindOrBuildElement("Fe");
    G4double density_SS = 8.06*CLHEP::g/CLHEP::cm3;
    G4int ncomponents_SS =6;
    G4double fractionmass;

    G4Material* StainlessSteel = 
            new G4Material("StainlessSteel", density_SS, ncomponents_SS);
    StainlessSteel->AddElement(elC, fractionmass=0.001);
    StainlessSteel->AddElement(elSi, fractionmass=0.007);
    StainlessSteel->AddElement(elCr, fractionmass=0.18);
    StainlessSteel->AddElement(elMn, fractionmass=0.01);
    StainlessSteel->AddElement(elFe, fractionmass=0.712);
    StainlessSteel->AddElement(elNi, fractionmass=0.09);

    // Visualization attributes
    G4VisAttributes* fBeamPipeVisAttribute =
            new G4VisAttributes(G4Colour(0.0,1.0,0.0));
    fBeamPipeVisAttribute->SetForceSolid(true);

    G4VisAttributes* fBeamPipeInsideVisAttribute =
            new G4VisAttributes(G4Colour(0.0,0.0,1.0));
    fBeamPipeInsideVisAttribute->SetForceSolid(false);

    // Variables
    G4double fMylarThickness = 10. * CLHEP::millimeter;

    //** BeamPipe **//
    G4Tubs* fBeamPipeEnvelopeSolid =
            new G4Tubs("BeamPipeEnvelope",
                       0.,
                       fBeamPipeRadius + fBeamPipeThickness,
                       length * 0.5 + fMylarThickness * 0.5 * 4.,
                       0*CLHEP::deg,
                       360*CLHEP::deg);
                   
    G4LogicalVolume* fBeamPipeEnvelopeLogic =
            new G4LogicalVolume(fBeamPipeEnvelopeSolid,
                                 Galactic,
                                 "BeamPipeEnvelope");


    G4Tubs* fBeamPipeSolid = new G4Tubs("BeamPipe",
                                        fBeamPipeRadius,
                                        fBeamPipeRadius + fBeamPipeThickness,
                                        length * 0.5,
                                        0*CLHEP::deg,
                                        360*CLHEP::deg);
    G4LogicalVolume* fBeamPipeLogic = new G4LogicalVolume(fBeamPipeSolid,
                                                          StainlessSteel,
                                                          "BeamPipe");
    fBeamPipeLogic->SetVisAttributes(fBeamPipeVisAttribute);
    new G4PVPlacement(0,
                      G4ThreeVector(),
                      fBeamPipeLogic,
                      "BeamPipe",
                      fBeamPipeEnvelopeLogic,
                      false,
                      0);


    G4Tubs* fBeamPipeInsideSolid = new G4Tubs("BeamPipeInside",
                                                0.,
                                                fBeamPipeRadius,
                                                length * 0.5,
                                                0*CLHEP::deg,
                                                360*CLHEP::deg);
    G4LogicalVolume* fBeamPipeInsideLogic =
            new G4LogicalVolume(fBeamPipeInsideSolid,
                            Vacuum,
                            "BeamPipeInside");
    fBeamPipeInsideLogic->SetVisAttributes(fBeamPipeInsideVisAttribute);
    new G4PVPlacement(0,
                      G4ThreeVector(),
                      fBeamPipeInsideLogic,
                      "BeamPipeInside",
                      fBeamPipeEnvelopeLogic,
                      false,
                      0);


    G4Tubs* fBeamPipeMylarSolid =
            new G4Tubs("BeamPipeMylar",
                       0.,
                       fBeamPipeRadius + fBeamPipeThickness,
                   fMylarThickness * 0.5,
                   0*CLHEP::deg,
                   360*CLHEP::deg);
                   
    G4LogicalVolume* fBeamPipeMylarLogic =
            new G4LogicalVolume(fBeamPipeMylarSolid,
                                Mylar,
                                "BeamPipeMylar");
    fBeamPipeMylarLogic->SetVisAttributes(fBeamPipeInsideVisAttribute);
    new G4PVPlacement(0,
                          G4ThreeVector(0.,
                                        0.,
                                        +(length + fMylarThickness * 2.) / 2.),
                          fBeamPipeMylarLogic,
                      "BeamPipeMylar",
                      fBeamPipeEnvelopeLogic,
                      false,
                      0);
                      
    new G4PVPlacement(0,
                          G4ThreeVector(0.,
                                        0.,
                                        -(length + fMylarThickness * 2.) / 2.),
                          fBeamPipeMylarLogic,
                      "BeamPipeMylar",
                      fBeamPipeEnvelopeLogic,
                      false,
                      1);

    return fBeamPipeEnvelopeLogic;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChDetectorConstruction::ConstructXtalTarget(){
    if(fXtalCurvatureRadius.x() != 0.){
        double fXtalAngleOut =
            fXtalAngle.y() + fXtalSize.z()/fXtalCurvatureRadius.x();
        fXtalSolid = new G4Tubs("Target",
                                fXtalCurvatureRadius.x() - fXtalSize.x()/2,
                                fXtalCurvatureRadius.x() + fXtalSize.x()/2,
                                fXtalSize.y()/2.,
                                fXtalAngle.y(),
                                fXtalAngleOut);
    }
    else{
        fXtalSolid = new G4Box("Target",
                               fXtalSize.x()/2.,
                               fXtalSize.y()/2.,
                               fXtalSize.z()/2.);
    }
    
    fXtalLogic = new G4LogicalVolume(fXtalSolid,fXtalMaterial,"Target");
    
    G4RotationMatrix* vRotationMatrix =
        new G4RotationMatrix; // Rotates X and Z axes only
    
    if(fXtalCurvatureRadius.x() != 0.){
        vRotationMatrix->rotateX(fXtalAngle.x()-CLHEP::pi*0.5);
        vRotationMatrix->rotateY(fXtalAngle.y());
        vRotationMatrix->rotateZ(fXtalAngle.z());
        fXtalPhysical =
            new G4PVPlacement(vRotationMatrix,
                              G4ThreeVector(-fXtalCurvatureRadius.x(),0.,0.),
                              fXtalLogic,"Target",
                              fWorldLogic,
                              false,
                              0);
    }
    else{
        vRotationMatrix->rotateX(fXtalAngle.x());
        vRotationMatrix->rotateY(fXtalAngle.y());
        vRotationMatrix->rotateZ(fXtalAngle.z());
        fXtalPhysical =
            new G4PVPlacement(vRotationMatrix,
                              G4ThreeVector(0.,0.,0.),
                              fXtalLogic,"Target",
                              fWorldLogic,
                              false,
                              0);
    }
    
    //----------------------------------------
    // Create XLogicalLattice
    //----------------------------------------
    XLogicalLattice* logicalLattice = new XLogicalLattice();
    double vScatteringConstant =
        3.67e-41*CLHEP::second*CLHEP::second*CLHEP::second;
    logicalLattice->SetScatteringConstant(vScatteringConstant);
    
    //----------------------------------------
    // Create XLogicalBase
    //----------------------------------------
    XLogicalAtomicLatticeDiamond *diamond_lattice =
        new XLogicalAtomicLatticeDiamond();
    G4Element* element = G4NistManager::
        Instance()->FindOrBuildElement(G4int(fXtalMaterial->GetZ()));
    XLogicalBase *base = new XLogicalBase(element,diamond_lattice);
    
    //----------------------------------------
    // Create XUnitCell
    //----------------------------------------
    XUnitCell* myCell = new XUnitCell();
    myCell->SetSize(fXtalCellSize);
    myCell->AddBase(base);
    
    //----------------------------------------
    // Create XPhysicalLattice
    //----------------------------------------
    XPhysicalLattice* physicalLattice =
        new XPhysicalLattice(fXtalPhysical, logicalLattice);
    physicalLattice->SetUnitCell(myCell);
    physicalLattice->SetMillerOrientation(G4int(fXtalMiller.x()),
                                          G4int(fXtalMiller.y()),
                                          G4int(fXtalMiller.z()));
    physicalLattice->SetLatticeOrientation(fXtalAngle.x(),
                                           fXtalAngle.y(),
                                           fXtalAngle.z());
    physicalLattice->SetThermalVibrationAmplitude(fXtalTVA);
    physicalLattice->SetCurvatureRadius(fXtalCurvatureRadius);
    
    //----------------------------------------
    // Register XPhysicalLattice
    //----------------------------------------
    if(XLatticeManager3::
        GetXLatticeManager()->GetXPhysicalLattice(fXtalPhysical)
        != physicalLattice){
            XLatticeManager3::
                GetXLatticeManager()->RegisterLattice(physicalLattice);
    }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.

void ExExChDetectorConstruction::SetXtalMaterial(const G4String& name){
    G4Material* vMaterial = G4Material::GetMaterial(name, false);
    
    if(!vMaterial){
        vMaterial = G4NistManager::Instance()->FindOrBuildMaterial(name);
    }
    
    if (vMaterial && vMaterial != fXtalMaterial) {
        G4cout << "DetectorConstructor::SetXtalMaterial() - New Xtal Material: "
            << vMaterial->GetName() << G4endl;
        fXtalMaterial = vMaterial;
        if(fXtalLogic){
            fXtalLogic->SetMaterial(vMaterial);
            G4RunManager::GetRunManager()->PhysicsHasBeenModified();
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String ExExChDetectorConstruction::GetXtalMaterial(){
    if(fXtalMaterial) {
        return fXtalMaterial->GetName();
    }
    return "";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChDetectorConstruction::SetXtalCurvatureRadius(G4ThreeVector cr){
    if(fXtalCurvatureRadius != cr) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fXtalCurvatureRadius = cr;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChDetectorConstruction::SetXtalSize(G4ThreeVector size) {
    if(fXtalSize != size) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fXtalSize = size;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChDetectorConstruction::SetXtalAngle(G4ThreeVector angle) {
    if(fXtalAngle != angle) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fXtalAngle = angle;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChDetectorConstruction::SetXtalCellSize(G4ThreeVector cellsize) {
    if(fXtalCellSize != cellsize) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fXtalCellSize = cellsize;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChDetectorConstruction::SetXtalMiller(G4ThreeVector miller) {
    if(fXtalMiller != miller) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fXtalMiller = miller;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChDetectorConstruction::SetXtalCellAngle(
                                            G4ThreeVector cellangle) {
    if(fXtalCellAngle != cellangle) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fXtalCellAngle = cellangle;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChDetectorConstruction::SetXtalThermalVibrationAmplitude(
                                            G4double thermvibr) {
    if(fXtalTVA != thermvibr) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fXtalTVA = thermvibr;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..
