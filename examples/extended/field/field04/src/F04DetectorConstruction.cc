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
/// \file field/field04/src/F04DetectorConstruction.cc
/// \brief Implementation of the F04DetectorConstruction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"

#include "F04DetectorConstruction.hh"
#include "F04DetectorMessenger.hh"

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "F04GlobalField.hh"

#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4RunManager.hh"
#include "G4StateManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "F04Materials.hh"

#include "G4RotationMatrix.hh"

#include "F04SimpleSolenoid.hh"
#include "F04FocusSolenoid.hh"

#include "G4AutoDelete.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04DetectorConstruction::F04DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fDetectorMessenger(0),
   fSolidWorld(0), fLogicWorld(0), fPhysiWorld(0),
   fSolidTarget(0), fLogicTarget(0), fPhysiTarget(0),
   fSolidDegrader(0), fLogicDegrader(0), fPhysiDegrader(0),
   fSolidCaptureMgnt(0), fLogicCaptureMgnt(0), fPhysiCaptureMgnt(0),
   fSolidTransferMgnt(0), fLogicTransferMgnt(0), fPhysiTransferMgnt(0),
   fWorldMaterial(0), fTargetMaterial(0), fDegraderMaterial(0)
{
  fWorldSizeZ = 50.*m;
  fWorldSizeR =  5.*m;

  fTargetRadius       =   0.4*cm;
  fTargetThickness    =  16.0*cm;

  SetTargetAngle(170);

  fDegraderRadius     =  30.0*cm;
  fDegraderThickness  =   0.1*cm;

  fCaptureMgntRadius  =  0.6*m;
  fCaptureMgntLength  =  4.0*m;

  SetCaptureMgntB1(2.5*tesla);
  SetCaptureMgntB2(5.0*tesla);

  fTransferMgntRadius =  0.3*m;
  fTransferMgntLength = 15.0*m;

  SetTransferMgntB(5.0*tesla);

  fDegraderPos = -fTransferMgntLength/2. + fDegraderThickness/2.;

  fDetectorMessenger = new F04DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04DetectorConstruction::~F04DetectorConstruction()
{
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* F04DetectorConstruction::Construct()
{

  if (fPhysiWorld) {
     G4GeometryManager::GetInstance()->OpenGeometry();
     G4PhysicalVolumeStore::GetInstance()->Clean();
     G4LogicalVolumeStore::GetInstance()->Clean();
     G4SolidStore::GetInstance()->Clean();
  }

  fMaterials = F04Materials::GetInstance();

  DefineMaterials();

  return ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorConstruction::DefineMaterials()
{
  //define materials for the experiment

  fVacuum = fMaterials->GetMaterial("G4_Galactic");

  fWorldMaterial    = fMaterials->GetMaterial("G4_AIR");
  fDegraderMaterial = fMaterials->GetMaterial("G4_Pb");
  fTargetMaterial   = fMaterials->GetMaterial("G4_W");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* F04DetectorConstruction::ConstructDetector()
{
  fSolidWorld = new G4Tubs("World",
               0.,GetWorldSizeR(),GetWorldSizeZ()/2.,0.,twopi);
 
  fLogicWorld = new G4LogicalVolume(fSolidWorld,
                                   GetWorldMaterial(),
                                   "World");

  fPhysiWorld = new G4PVPlacement(0,
                                 G4ThreeVector(),
                                 "World",
                                 fLogicWorld,
                                 0,
                                 false,
                                 0);

  // Capture Magnet

  fSolidCaptureMgnt = new G4Tubs("CaptureMgnt",
                     0.,GetCaptureMgntRadius(),
                     GetCaptureMgntLength()/2.,0.,twopi);

  fLogicCaptureMgnt = new G4LogicalVolume(fSolidCaptureMgnt,
                                         fVacuum,
                                         "CaptureMgnt");

  fCaptureMgntCenter = G4ThreeVector();

  fPhysiCaptureMgnt = new G4PVPlacement(0,
                                       fCaptureMgntCenter,
                                       "CaptureMgnt",
                                       fLogicCaptureMgnt,
                                       fPhysiWorld,
                                       false,
                                       0);

  // Transfer Magnet

  fSolidTransferMgnt = new G4Tubs("TransferMgnt",
                      0.,GetTransferMgntRadius(),
                      GetTransferMgntLength()/2.,0.,twopi);

  fLogicTransferMgnt = new G4LogicalVolume(fSolidTransferMgnt,
                                          fVacuum,
                                          "TransferMgnt");

  G4double z  = GetCaptureMgntLength()/2. + GetTransferMgntLength()/2.
                              + GetTransferMgntPos();
  G4double x =  GetTransferMgntPos()/2.;

  fTransferMgntCenter = G4ThreeVector(x,0.,z);

  G4RotationMatrix* g4rot = new G4RotationMatrix();
  *g4rot = StringToRotationMatrix("Y30,X10");
  *g4rot = g4rot->inverse();
  if (*g4rot == G4RotationMatrix()) g4rot = NULL;

  fPhysiTransferMgnt = new G4PVPlacement(g4rot,
                                        fTransferMgntCenter,
                                        "TransferMgnt",
                                        fLogicTransferMgnt,
                                        fPhysiWorld,
                                        false,
                                        0);

  // Test Plane

  G4Tubs* solidTestPlane = new G4Tubs("TestPlane",
                                      0.,GetTransferMgntRadius(),
                                      1.*mm,0.,twopi);

  G4LogicalVolume* logicTestPlane = new G4LogicalVolume(solidTestPlane,
                                                        fVacuum,
                                                        "TestPlane");

  z = GetTransferMgntLength()/2. - 1.*mm;

  G4ThreeVector testPlaneCenter = G4ThreeVector(0.,0.,z);

  new G4PVPlacement(0,
                    testPlaneCenter,
                    "TestPlane",
                    logicTestPlane,
                    fPhysiTransferMgnt,
                    false,
                    0);

  // Target

  if (GetTargetThickness() > 0.)
  {
      fSolidTarget = new G4Tubs("Target",
                    0.,GetTargetRadius(),
                    GetTargetThickness()/2.,0.,twopi);

      fLogicTarget = new G4LogicalVolume(fSolidTarget,
                                        GetTargetMaterial(),
                                        "Target");

      G4int i =  GetTargetAngle();

      G4String angle = std::to_string(i);
      G4StrUtil::strip(angle);
      angle = "Y" + angle;

      g4rot = new G4RotationMatrix();
      *g4rot = StringToRotationMatrix(angle);
      *g4rot = g4rot->inverse();
      if (*g4rot == G4RotationMatrix()) g4rot = NULL;

      G4ThreeVector targetCenter(0.,0.,GetTargetPos());

      fPhysiTarget = new G4PVPlacement(g4rot,
                                      targetCenter,
                                      "Target",
                                      fLogicTarget,
                                      fPhysiCaptureMgnt,
                                      false,
                                      0);
  }

  // Degrader

  if (GetDegraderThickness() > 0.)
  {
      fSolidDegrader = new G4Tubs("Degrader",
                      0., GetDegraderRadius(),
                      GetDegraderThickness()/2., 0.,twopi);
 
      fLogicDegrader = new G4LogicalVolume(fSolidDegrader,
                                                GetDegraderMaterial(),
                                                "Degrader");

      G4ThreeVector degraderCenter = G4ThreeVector(0.,0.,GetDegraderPos());

      fPhysiDegrader = new G4PVPlacement(0,
                                        degraderCenter,
                                        "Degrader",
                                        fLogicDegrader,
                                        fPhysiTransferMgnt,
                                        false,
                                        0);
  }

  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorConstruction::SetWorldMaterial(const G4String materialChoice)
{
  G4Material* pttoMaterial =
      G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial != fWorldMaterial) {
     if ( pttoMaterial ) {
        fWorldMaterial = pttoMaterial;
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
     } else {
        G4cout << "\n--> WARNING from SetWorldMaterial : "
               << materialChoice << " not found" << G4endl;
     }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorConstruction::SetTargetMaterial(const G4String materialChoice)
{
  G4Material* pttoMaterial =
      G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial != fTargetMaterial) {
     if ( pttoMaterial ) {
        fTargetMaterial = pttoMaterial;
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
     } else {
        G4cout << "\n-->  WARNING from SetTargetMaterial : "
               << materialChoice << " not found" << G4endl;
     }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorConstruction::SetDegraderMaterial(const G4String materialChoice)

{
  G4Material* pttoMaterial =
      G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial != fDegraderMaterial) {
     if ( pttoMaterial ) {
        fDegraderMaterial = pttoMaterial;
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
     } else {
        G4cout << "\n--> WARNING from SetDegraderMaterial : "
               << materialChoice << " not found" << G4endl;
     }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorConstruction::SetWorldSizeZ(G4double val)
{
  fWorldSizeZ = val;

  if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorConstruction::SetWorldSizeR(G4double val)
{
  fWorldSizeR = val;
  if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorConstruction::SetCaptureMgntRadius(G4double val)
{
  fCaptureMgntRadius = val;
  if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorConstruction::SetCaptureMgntLength(G4double val)
{
  fCaptureMgntLength = val;
  if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorConstruction::SetCaptureMgntB1(G4double val)
{
  fCaptureMgntB1 = val;
  if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorConstruction::SetCaptureMgntB2(G4double val)
{
  fCaptureMgntB2 = val;
  if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorConstruction::SetTransferMgntRadius(G4double val)
{
  fTransferMgntRadius = val;
  if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorConstruction::SetTransferMgntLength(G4double val)
{
  fTransferMgntLength = val;
  if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorConstruction::SetTransferMgntB(G4double val)
{
  fTransferMgntB = val;
  if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorConstruction::SetTransferMgntPos(G4double val)
{
  fTransferMgntPos = val;
  if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorConstruction::SetTargetRadius(G4double val)
{
  fTargetRadius = val;
  if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorConstruction::SetTargetThickness(G4double val)
{
  fTargetThickness = val;
  if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorConstruction::SetTargetPos(G4double val)
{
  fTargetPos = val;
  if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorConstruction::SetTargetAngle(G4int val)
{
  fTargetAngle = val;
  if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorConstruction::SetDegraderRadius(G4double val)
{
  fDegraderRadius = val;
  if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorConstruction::SetDegraderThickness(G4double val)
{
  fDegraderThickness = val;
  if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorConstruction::SetDegraderPos(G4double val)
{
  fDegraderPos = val;
  if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorConstruction::ConstructSDandField()
{
  if (!fFieldSetUp.Get()) {
     F04GlobalField* field = F04GlobalField::GetObject(this);
     G4AutoDelete::Register(field);  // Kernel will delete the F04GlobalField
     fFieldSetUp.Put(field);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4RotationMatrix
            F04DetectorConstruction::StringToRotationMatrix(G4String rotation)
{
  // We apply successive rotations OF THE OBJECT around the FIXED
  // axes of the parent's local coordinates; rotations are applied
  // left-to-right (rotation="r1,r2,r3" => r1 then r2 then r3).

  G4RotationMatrix rot;

  unsigned int place = 0;

  while (place < rotation.size()) {

        G4double angle;
        char* p(0);
        G4String current=rotation.substr(place+1);
        angle = strtod(current.c_str(),&p) * deg;

        if (!p || (*p != ',' && *p != '\0')) {
          G4cerr << "Invalid rotation specification: " <<
                                                  rotation.c_str() << G4endl;

           return rot;
        }

        G4RotationMatrix thisRotation;

        switch(rotation.substr(place,1).c_str()[0]) {
              case 'X': case 'x':
                thisRotation = G4RotationMatrix(CLHEP::HepRotationX(angle));
                break;
              case 'Y': case 'y':
                thisRotation = G4RotationMatrix(CLHEP::HepRotationY(angle));
                break;
              case 'Z': case 'z':
                thisRotation = G4RotationMatrix(CLHEP::HepRotationZ(angle));
                break;
              default:
                G4cerr << " Invalid rotation specification: "
                       << rotation << G4endl;
                return rot;
        }

       rot = thisRotation * rot;
       place = rotation.find(',',place);
       if (place > rotation.size()) break;
       ++place;
  }

  return rot;

}
