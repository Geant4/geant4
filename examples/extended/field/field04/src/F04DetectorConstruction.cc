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
//

#include "G4ios.hh"
#include "globals.hh"

#include "F04DetectorConstruction.hh"
#include "F04DetectorMessenger.hh"

#include "F04GlobalField.hh"

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"
#include "G4TransportationManager.hh"

#include "G4GeometryManager.hh"

#include "G4SolidStore.hh"
#include "G4RegionStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4RunManager.hh"

#include "F04DetectorConstruction.hh"
#include "F04DetectorMessenger.hh"
#include "F04Materials.hh"

#include "G4RotationMatrix.hh"

#include "F04SimpleSolenoid.hh"
#include "F04FocusSolenoid.hh"

F04DetectorConstruction::F04DetectorConstruction()
 : solidWorld(0), logicWorld(0), physiWorld(0),
   solidTarget(0), logicTarget(0), physiTarget(0),
   solidDegrader(0),logicDegrader(0), physiDegrader(0),
   solidCaptureMgnt(0), logicCaptureMgnt(0), physiCaptureMgnt(0),
   solidTransferMgnt(0), logicTransferMgnt(0), physiTransferMgnt(0),
   WorldMaterial(0), TargetMaterial(0), DegraderMaterial(0),
   focusSolenoid(0), simpleSolenoid(0)
{
  WorldSizeZ = 50.*m;
  WorldSizeR =  5.*m;

  TargetRadius       =   0.4*cm;
  TargetThickness    =  16.0*cm;

  SetTargetAngle(170);

  DegraderRadius     =  30.0*cm;
  DegraderThickness  =   0.1*cm;

  CaptureMgntRadius  =  0.6*m;
  CaptureMgntLength  =  4.0*m;

  SetCaptureMgntB1(2.5*tesla);
  SetCaptureMgntB2(5.0*tesla);

  TransferMgntRadius =  0.3*m;
  TransferMgntLength = 15.0*m;

  SetTransferMgntB(5.0*tesla);

  DegraderPos = -TransferMgntLength/2. + DegraderThickness/2.;

  // ensure the global field is initialized
  (void)F04GlobalField::getObject();

  detectorMessenger = new F04DetectorMessenger(this);
}

F04DetectorConstruction::~F04DetectorConstruction()
{ 
  delete detectorMessenger;
}

G4VPhysicalVolume* F04DetectorConstruction::Construct()
{
  materials = F04Materials::GetInstance();

  DefineMaterials();

  return ConstructDetector();
}

void F04DetectorConstruction::DefineMaterials()
{
  //define materials for the experiment

  Vacuum = materials->GetMaterial("G4_Galactic");

  WorldMaterial    = materials->GetMaterial("G4_AIR");
  DegraderMaterial = materials->GetMaterial("G4_Pb");
  TargetMaterial   = materials->GetMaterial("G4_W");
}

G4VPhysicalVolume* F04DetectorConstruction::ConstructDetector()
{
  solidWorld = new G4Tubs("World",
               0.,GetWorldSizeR(),GetWorldSizeZ()/2.,0.,twopi);
                         
  logicWorld = new G4LogicalVolume(solidWorld,
                                   GetWorldMaterial(),
                                   "World");
                                   
  physiWorld = new G4PVPlacement(0,
  				 G4ThreeVector(),
                                 "World",
                                 logicWorld,
                                 0,
                                 false,
                                 0);

  // Capture Magnet

  solidCaptureMgnt = new G4Tubs("CaptureMgnt",
                     0.,GetCaptureMgntRadius(),
                     GetCaptureMgntLength()/2.,0.,twopi);

  logicCaptureMgnt = new G4LogicalVolume(solidCaptureMgnt,
                                         Vacuum,
                                         "CaptureMgnt");

  G4ThreeVector CaptureMgntCenter = G4ThreeVector();

  physiCaptureMgnt = new G4PVPlacement(0,
                                       CaptureMgntCenter,
                                       "CaptureMgnt",
                                       logicCaptureMgnt,
                                       physiWorld,
                                       false,
                                       0);

  // Transfer Magnet

  solidTransferMgnt = new G4Tubs("TransferMgnt",
                      0.,GetTransferMgntRadius(),
                      GetTransferMgntLength()/2.,0.,twopi);

  logicTransferMgnt = new G4LogicalVolume(solidTransferMgnt,
                                          Vacuum,
                                          "TransferMgnt");

  G4double z  = GetCaptureMgntLength()/2. + GetTransferMgntLength()/2.
                              + GetTransferMgntPos();
  G4double x =  GetTransferMgntPos()/2.;

  G4ThreeVector TransferMgntCenter = G4ThreeVector(x,0.,z);

  G4RotationMatrix* g4rot = new G4RotationMatrix();
  *g4rot = stringToRotationMatrix("Y30,X10");
  *g4rot = g4rot->inverse();
  if (*g4rot == G4RotationMatrix()) g4rot = NULL;

  physiTransferMgnt = new G4PVPlacement(g4rot,
                                        TransferMgntCenter,
                                        "TransferMgnt",
                                        logicTransferMgnt,
                                        physiWorld,
                                        false,
                                        0);

  // Test Plane

  G4Tubs* solidTestPlane = new G4Tubs("TestPlane",
                                      0.,GetTransferMgntRadius(),
                                      1.*mm,0.,twopi);

  G4LogicalVolume* logicTestPlane = new G4LogicalVolume(solidTestPlane,
                                                        Vacuum,
                                                        "TestPlane");


  z = GetTransferMgntLength()/2. - 1.*mm;

  G4ThreeVector TestPlaneCenter = G4ThreeVector(0.,0.,z);

  new G4PVPlacement(0,
                    TestPlaneCenter,
                    "TestPlane",
                    logicTestPlane,
                    physiTransferMgnt,
                    false,
                    0);

  // Target

  if (GetTargetThickness() > 0.)
  {
      solidTarget = new G4Tubs("Target", 
                    0.,GetTargetRadius(),
                    GetTargetThickness()/2.,0.,twopi);

      logicTarget = new G4LogicalVolume(solidTarget,
                                        GetTargetMaterial(),
                                        "Target");

      G4int i =  GetTargetAngle();

      char c[4];
      sprintf(c,"%d",i);
      G4String Angle = c;
      Angle = Angle.strip(G4String::both,' ');
      Angle = "Y" + Angle;

      g4rot = new G4RotationMatrix();
      *g4rot = stringToRotationMatrix(Angle);
      *g4rot = g4rot->inverse();
      if (*g4rot == G4RotationMatrix()) g4rot = NULL;

      G4ThreeVector TargetCenter(0.,0.,GetTargetPos());

      physiTarget = new G4PVPlacement(g4rot,
                                      TargetCenter,
                                      "Target",
                                      logicTarget,
                                      physiCaptureMgnt,
                                      false,
                                      0);
  }

  // Degrader

  if (GetDegraderThickness() > 0.) 
  { 
      solidDegrader = new G4Tubs("Degrader",
                      0., GetDegraderRadius(),
                      GetDegraderThickness()/2., 0.,twopi); 
                          
      logicDegrader = new G4LogicalVolume(solidDegrader,    
      			                  GetDegraderMaterial(), 
      			                  "Degrader");     

      G4ThreeVector DegraderCenter = G4ThreeVector(0.,0.,GetDegraderPos());

      physiDegrader = new G4PVPlacement(0,		   
      		                        DegraderCenter,
                                        "Degrader",        
                                        logicDegrader,     
                                        physiTransferMgnt,       
                                        false,             
                                        0);
  }

  G4double l = 0.0;
  G4double B1 = GetCaptureMgntB1();
  G4double B2 = GetCaptureMgntB2();

  if (focusSolenoid) delete focusSolenoid;
  focusSolenoid = new F04FocusSolenoid(B1, B2, l,
                                    logicCaptureMgnt,CaptureMgntCenter);
  focusSolenoid -> SetHalf(true);

           l = 0.0;
  G4double B = GetTransferMgntB();

  if (simpleSolenoid) delete simpleSolenoid;
  simpleSolenoid = new F04SimpleSolenoid(B, l,
                                      logicTransferMgnt,TransferMgntCenter);

  simpleSolenoid->setColor("1,0,1");
  simpleSolenoid->setColor("0,1,1");
  simpleSolenoid->setMaxStep(1.5*mm);
  simpleSolenoid->setMaxStep(2.5*mm);

  return physiWorld;
}

void F04DetectorConstruction::SetWorldMaterial(const G4String materialChoice)
{
  G4Material* pttoMaterial =
      G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial != WorldMaterial) {
     if ( pttoMaterial ) {
        WorldMaterial = pttoMaterial;
     } else {
        G4cout << "\n--> WARNING from SetWorldMaterial : "
               << materialChoice << " not found" << G4endl;
     }
  }
}

void F04DetectorConstruction::SetTargetMaterial(const G4String materialChoice)
{
  G4Material* pttoMaterial =
      G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial != TargetMaterial) {
     if ( pttoMaterial ) {
        TargetMaterial = pttoMaterial;
     } else {
        G4cout << "\n-->  WARNING from SetTargetMaterial : "
               << materialChoice << " not found" << G4endl;
     }
  }
}

void F04DetectorConstruction::SetDegraderMaterial(const G4String materialChoice)

{
  G4Material* pttoMaterial =
      G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial != DegraderMaterial) {
     if ( pttoMaterial ) {
        DegraderMaterial = pttoMaterial;
     } else {
        G4cout << "\n--> WARNING from SetDegraderMaterial : "
               << materialChoice << " not found" << G4endl;
     }
  }
}

void F04DetectorConstruction::SetWorldSizeZ(G4double val)
{
  WorldSizeZ = val;
}

void F04DetectorConstruction::SetWorldSizeR(G4double val)
{
  WorldSizeR = val;
}

void F04DetectorConstruction::SetCaptureMgntRadius(G4double val)
{
  CaptureMgntRadius = val;
}

void F04DetectorConstruction::SetCaptureMgntLength(G4double val)
{
  CaptureMgntLength = val;
}

void F04DetectorConstruction::SetCaptureMgntB1(G4double val)
{
  CaptureMgntB1 = val;
}

void F04DetectorConstruction::SetCaptureMgntB2(G4double val)
{
  CaptureMgntB2 = val;
}

void F04DetectorConstruction::SetTransferMgntRadius(G4double val)
{
  TransferMgntRadius = val;
}

void F04DetectorConstruction::SetTransferMgntLength(G4double val)
{
  TransferMgntLength = val;
}

void F04DetectorConstruction::SetTransferMgntB(G4double val)
{
  TransferMgntB = val;
}

void F04DetectorConstruction::SetTransferMgntPos(G4double val)
{
  TransferMgntPos = val;
}

void F04DetectorConstruction::SetTargetRadius(G4double val)
{
  TargetRadius = val;
}

void F04DetectorConstruction::SetTargetThickness(G4double val)
{
  TargetThickness = val;
}

void F04DetectorConstruction::SetTargetPos(G4double val)
{
  TargetPos = val;
}

void F04DetectorConstruction::SetTargetAngle(G4int val)
{
  TargetAngle = val;
}

void F04DetectorConstruction::SetDegraderRadius(G4double val)
{
  DegraderRadius = val;
}

void F04DetectorConstruction::SetDegraderThickness(G4double val)
{
  DegraderThickness = val;
}  

void F04DetectorConstruction::SetDegraderPos(G4double val)
{
  DegraderPos = val;
}  

void F04DetectorConstruction::UpdateGeometry()
{
  if (!physiWorld) return;

  // clean-up previous geometry
  G4GeometryManager::GetInstance()->OpenGeometry();

  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  //define new one
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetector());

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();

  G4RegionStore::GetInstance()->UpdateMaterialList(physiWorld);

}

G4RotationMatrix 
            F04DetectorConstruction::stringToRotationMatrix(G4String rotation)
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
