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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4ReflectionFactory.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fMessenger(0),
   fMethod(kWithDirectMatrix),
   fWorldVolume(0),
   fTrdVolume(0)
   
{
  fMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Materials
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* material = nist->FindOrBuildMaterial("G4_AIR");
  
  // Clean old geometry, if any
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
        
  // World
  //  
  G4double rmin = 0.;
  G4double rmax = 5*cm;
  G4double hz   = 5*cm;
  G4double phiMin = 0.;
  G4double deltaPhi = 360*degree;
        
  G4Tubs* solidWorld 
    = new G4Tubs("World",                            //name
                  rmin, rmax, hz, phiMin, deltaPhi); //size

  fWorldVolume 
    = new G4LogicalVolume(solidWorld,                //solid
                          material,                  //material
                          "World");                  //name

  G4VPhysicalVolume* physiWorld
    = new G4PVPlacement(0,                     //no rotation
                        G4ThreeVector(),       //at (0,0,0)
                        fWorldVolume,          //logical volume
                        "World",               //name
                        0,                     //mother volume
                        false,                 //no boolean operation
                        0);                    //copy number

  // Trd volume 
  //  
  G4double dX1 = 1*cm;
  G4double dX2 = 1*cm;
  G4double dY1 = 1*cm; 
  G4double dY2 = 2*cm;
  G4double dZ  = 3*cm; 

  G4Trd* solidTrd 
    = new G4Trd("trd",                               //name
                 dX1/2, dX2/2, dY1/2, dY2/2, dZ/2);  //size

  fTrdVolume
    = new G4LogicalVolume(solidTrd,                  //solid
                          material,                  //material
                          "trd");                    //name
                                

  // Place Volume1 and Volume2 according to selected methods
  // 
  switch ( fMethod ) {
    case kWithDirectMatrix:   PlaceWithDirectMatrix(); break;
    case kWithInverseMatrix:  PlaceWithInverseMatrix(); break;
    case kWithAxialRotations: PlaceWithAxialRotations(); break;
    case kWithEulerAngles:    PlaceWithEulerAngles(); break;
    case kWithReflections:    PlaceWithReflections(); break;
    default: ;;
  }

  // Return the root volume
  //
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PlaceWithDirectMatrix()
{
  G4double og = 3*cm; 

  // 1st position 
  //
  G4double phi = 30*deg;
  // u, v, w are the daughter axes, projected on the mother frame
  G4ThreeVector u = G4ThreeVector(0, 0, -1);
  G4ThreeVector v = G4ThreeVector(-std::sin(phi), std::cos(phi),0.);
  G4ThreeVector w = G4ThreeVector( std::cos(phi), std::sin(phi),0.);
  G4RotationMatrix rotm1  = G4RotationMatrix(u, v, w);
  G4cout << "\n --> phi = " << phi/deg << " deg;  direct rotation matrix : ";
  rotm1.print(G4cout);     
  G4ThreeVector position1 = og*w;
  G4Transform3D transform1 = G4Transform3D(rotm1,position1);

  new G4PVPlacement(transform1,         //position, rotation        
                    fTrdVolume,         //logical volume
                    "Trd",              //name
                    fWorldVolume,       //mother volume
                    false,              //no boolean operation
                    1);                 //copy number                       

  // 2nd position 
  //
  phi = phi + 90*deg;
  v = G4ThreeVector(-std::sin(phi), std::cos(phi),0.);
  w = G4ThreeVector( std::cos(phi), std::sin(phi),0.);
  G4RotationMatrix rotm2  = G4RotationMatrix(u, v, w);
  G4ThreeVector position2 = og*w;
  G4Transform3D transform2 = G4Transform3D(rotm2,position2);
  new G4PVPlacement(transform2,         //position, rotation        
                    fTrdVolume,         //logical volume
                    "Trd",              //name
                    fWorldVolume,       //mother volume
                    false,              //no boolean operation
                    2);                 //copy number                       
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PlaceWithInverseMatrix()
{
  G4double og = 3*cm; 

  // 1st position
  //
  G4double phi = 30*deg;
  // u, v, w are the daughter axes, projected on the mother frame     
  G4ThreeVector u = G4ThreeVector(0, 0, -1);
  G4ThreeVector v = G4ThreeVector(-std::sin(phi), std::cos(phi),0.);
  G4ThreeVector w = G4ThreeVector( std::cos(phi), std::sin(phi),0.);
  G4RotationMatrix rotm1 = G4RotationMatrix(u, v, w); 
  G4RotationMatrix* rotm1Inv = new G4RotationMatrix(rotm1.inverse());
  G4cout << "\n --> phi = " << phi/deg << " deg;  inverse rotation matrix : ";
  rotm1Inv->print(G4cout);   
  G4ThreeVector position1 = og*w;

  new G4PVPlacement(rotm1Inv,
                    position1,
                    fTrdVolume,         //logical volume
                    "Trd",              //name
                    fWorldVolume,       //mother volume
                    false,              //no boolean operation
                    1);                 //copy number                       

  // 2nd position
  //
  phi = phi + 90*deg;
  v = G4ThreeVector(-std::sin(phi), std::cos(phi),0.);
  w = G4ThreeVector( std::cos(phi), std::sin(phi),0.);
  G4RotationMatrix rotm2  = G4RotationMatrix(u, v, w);
  G4RotationMatrix* rotm2Inv  = new G4RotationMatrix(rotm2.inverse());
  G4ThreeVector position2 = og*w;

  new G4PVPlacement(rotm2Inv,           //rotation
                    position2,          //position        
                    fTrdVolume,         //logical volume
                    "Trd",              //name
                    fWorldVolume,       //mother volume
                    false,              //no boolean operation
                    2);                 //copy number
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PlaceWithAxialRotations()
{
  G4double og = 3*cm; 

  // 1st position (with first G4PVPlacement constructor)
  //
  G4double phi = 30*deg, theta = 90*deg;
  G4ThreeVector rotAxis =
               G4ThreeVector(std::sin(theta-pi/2), 0., std::cos(theta-pi/2));
  G4RotationMatrix rotm1 = G4RotationMatrix();
  rotm1.rotateY(theta);
  rotm1.rotate (phi, rotAxis);  
  G4cout << "\n --> direct rotation matrix : "
         << " theta = " << theta/deg << " deg;"
         << " phi  = "  <<  phi/deg << " deg;";
  rotm1.print(G4cout);        
  G4ThreeVector w = G4ThreeVector( std::sin(theta)*std::cos(phi),
                    std::sin(phi), std::cos(theta)*std::cos(phi));
  G4ThreeVector position1 = og*w;
  G4Transform3D transform1(rotm1,position1);

  new G4PVPlacement(transform1,          //rotation,position
                    fTrdVolume,          //logical volume
                    "Trd",               //name
                    fWorldVolume,        //mother volume
                    false,               //no boolean operation
                    1);                  //copy number

  // 2nd position (with second G4PVPlacement constructor)
  //
  phi = phi + 90*deg;
  //rotm2Inv could be calculated with rotm2.inverse()
  //but also by the following :
  G4RotationMatrix* rotm2Inv  = new G4RotationMatrix();
  rotm2Inv->rotate (-phi, rotAxis);  
  rotm2Inv->rotateY(-theta);
  w = G4ThreeVector( std::sin(theta)*std::cos(phi),
                     std::sin(phi), std::cos(theta)*std::cos(phi));
  G4ThreeVector position2 = og*w;
  
  new G4PVPlacement(rotm2Inv,           //rotation
                    position2,          //position
                    fTrdVolume,         //logical volume
                    "Trd",              //name
                    fWorldVolume,       //mother volume
                    false,              //no boolean operation
                    2);                 //copy number
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PlaceWithEulerAngles()
{
  //definitions : mother frame = {x,y,z} ; daughter frame = {u,v,w}
  // n = node line = intercept of xy and uv planes
  // phi_euler   = (x,n) : precession
  // theta_euler = (z,w) : nutation
  // psi_euler   = (n,u) : proper rotation
  
  G4double og = 3*cm;  

  // 1st position (with first G4PVPlacement constructor)
  //
  G4double phi = 30*deg;   
  G4double phi_euler   =  phi + pi/2; 
  G4double theta_euler =  90*deg;
  G4double psi_euler   = -90*deg;
  //attention : clhep Euler constructor build inverse matrix !
  G4RotationMatrix rotm1Inv = G4RotationMatrix(phi_euler,theta_euler,psi_euler);
  G4RotationMatrix rotm1 = rotm1Inv.inverse();
  //remark : could be built as rotm1 = G4RotationMatrix(-psi, -theta, -phi)
  G4cout << "\n --> phi = " << phi/deg << " deg;  direct rotation matrix : ";
  rotm1.print(G4cout);           
  G4ThreeVector w = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);
  G4ThreeVector position1 = og*w;
  G4Transform3D transform1 = G4Transform3D(rotm1,position1);

  new G4PVPlacement(transform1,         //position, rotation        
                    fTrdVolume,         //logical volume
                    "Trd",              //name
                    fWorldVolume,       //mother volume
                    false,              //no boolean operation
                    1);                 //copy number                       

  // 2nd position (with second G4PVPlacement constructor)
  //
  phi = phi + 90*deg;
  
  phi_euler   =  phi + pi/2; 
  G4RotationMatrix* rotm2Inv  
                = new G4RotationMatrix(phi_euler,theta_euler,psi_euler);  
  w = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);
  G4ThreeVector position2 = og*w;
                           
  new G4PVPlacement(rotm2Inv,           //rotation
                    position2,          //position                     
                    fTrdVolume,         //logical volume
                    "Trd",              //name
                    fWorldVolume,       //mother volume
                    false,              //no boolean operation
                    2);                 //copy number                       
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PlaceWithReflections()
{
/// Placement with reflections.
/// In order to better show the reflection symmetry we do not apply
/// the rotation along Y axis.  

  G4double og = 3*cm; 

  // Place first two positionz in z = - 3cm
  //

  // 1st position
  G4double phi = 30*deg;
  G4RotationMatrix rotm1;
  //rotm1.rotateY(90*deg); 
  rotm1.rotateZ(phi); 
  G4ThreeVector uz = G4ThreeVector(std::cos(phi), std::sin(phi), 0);    
  G4ThreeVector position = og*uz;
  G4Transform3D transform1(rotm1,position);
  G4Transform3D translateZ = HepGeom::Translate3D(0, 0, -3.*cm);

  new G4PVPlacement(translateZ * transform1, //rotation,position
                    fTrdVolume,              //logical volume
                    "Trd",                   //name
                    fWorldVolume,            //mother volume
                    false,                   //no boolean operation
                    1);                      //copy number

  // 2nd position
  phi = phi + pi/2 ;
  G4RotationMatrix rotm2;
  //rotm2.rotateY(90*deg); 
  rotm2.rotateZ(phi); 
  uz = G4ThreeVector(std::cos(phi), std::sin(phi), 0.);    
  position = og*uz;
  G4Transform3D transform2 = G4Transform3D(rotm2, position);
  
  new G4PVPlacement(translateZ * transform2, //rotation, position
                    fTrdVolume,              //logical volume
                    "Trd",                   //name
                    fWorldVolume,            //mother volume
                    false,                   //no boolean operation
                    2);                      //copy number


  // Place next two positionz in z = + 3cm with reflection 
  //

  // 3rd position
  translateZ = HepGeom::Translate3D(0, 0, +3.*cm);
  G4Transform3D reflect3D = HepGeom::ReflectZ3D();

  G4ReflectionFactory::Instance()
    ->Place(translateZ * transform1 * reflect3D, //rotation,position
            "Trd",                               //name
            fTrdVolume,                          //logical volume
            fWorldVolume,                        //mother volume
            false,                               //no boolean operation
            3);                                  //copy number

  // 4rd position
  G4ReflectionFactory::Instance()
    ->Place( translateZ * transform2 * reflect3D,//rotation,position
            "Trd",                               //name
            fTrdVolume,                          //logical volume
            fWorldVolume,                        //mother volume
            false,                               //no boolean operation
            4);                                  //copy number
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"

void DetectorConstruction::SetMethod(EMethod method) { 
  fMethod = method;
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());  
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
