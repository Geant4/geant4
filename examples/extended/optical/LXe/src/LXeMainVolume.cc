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
/// \file optical/LXe/src/LXeMainVolume.cc
/// \brief Implementation of the LXeMainVolume class
//
//
#include "globals.hh"

#include "LXeMainVolume.hh"

#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeMainVolume::LXeMainVolume(G4RotationMatrix *pRot,
                             const G4ThreeVector &tlate,
                             G4LogicalVolume *pMotherLogical,
                             G4bool pMany,
                             G4int pCopyNo,
                             LXeDetectorConstruction* c)
  //Pass info to the G4PVPlacement constructor
  :G4PVPlacement(pRot,tlate,
                 //Temp logical volume must be created here
                 new G4LogicalVolume(new G4Box("temp",1,1,1),
                                     G4Material::GetMaterial("Vacuum"),
                                     "temp",0,0,0),
                 "housing",pMotherLogical,pMany,pCopyNo),fConstructor(c)
{
  CopyValues();

  G4double housing_x=fScint_x+2.*fD_mtl;
  G4double housing_y=fScint_y+2.*fD_mtl;
  G4double housing_z=fScint_z+2.*fD_mtl;
 
  //*************************** housing and scintillator
  fScint_box = new G4Box("scint_box",fScint_x/2.,fScint_y/2.,fScint_z/2.);
  fHousing_box = new G4Box("housing_box",housing_x/2.,housing_y/2.,
                           housing_z/2.);
 
  fScint_log = new G4LogicalVolume(fScint_box,G4Material::GetMaterial("LXe"),
                                   "scint_log",0,0,0);
  fHousing_log = new G4LogicalVolume(fHousing_box,
                                     G4Material::GetMaterial("Al"),
                                     "housing_log",0,0,0);
 
  new G4PVPlacement(0,G4ThreeVector(),fScint_log,"scintillator",
                                 fHousing_log,false,0);
 
  //*************** Miscellaneous sphere to demonstrate skin surfaces
  fSphere = new G4Sphere("sphere",0.*mm,2.*cm,0.*deg,360.*deg,0.*deg,360.*deg);
  fSphere_log = new G4LogicalVolume(fSphere,G4Material::GetMaterial("Al"),
                                    "sphere_log");
  if(fSphereOn)
    new G4PVPlacement(0,G4ThreeVector(5.*cm,5.*cm,5.*cm),
                                      fSphere_log,"sphere",fScint_log,false,0);
 
  //****************** Build PMTs
  G4double innerRadius_pmt = 0.*cm;
  G4double height_pmt = fD_mtl/2.;
  G4double startAngle_pmt = 0.*deg;
  G4double spanningAngle_pmt = 360.*deg;
 
  fPmt = new G4Tubs("pmt_tube",innerRadius_pmt,fOuterRadius_pmt,
                    height_pmt,startAngle_pmt,spanningAngle_pmt);
 
  //the "photocathode" is a metal slab at the back of the glass that
  //is only a very rough approximation of the real thing since it only
  //absorbs or detects the photons based on the efficiency set below
  fPhotocath = new G4Tubs("photocath_tube",innerRadius_pmt,fOuterRadius_pmt,
                          height_pmt/2,startAngle_pmt,spanningAngle_pmt);
 
  fPmt_log = new G4LogicalVolume(fPmt,G4Material::GetMaterial("Glass"),
                                 "pmt_log");
  fPhotocath_log = new G4LogicalVolume(fPhotocath,
                                       G4Material::GetMaterial("Al"),
                                       "photocath_log");
 
  new G4PVPlacement(0,G4ThreeVector(0,0,-height_pmt/2),
                                    fPhotocath_log,"photocath",
                                    fPmt_log,false,0);
 
  //***********Arrange pmts around the outside of housing**********

  G4double dx = fScint_x/fNx;
  G4double dy = fScint_y/fNy;
  G4double dz = fScint_z/fNz;
 
  G4double x,y,z;
  G4double xmin = -fScint_x/2. - dx/2.;
  G4double ymin = -fScint_y/2. - dy/2.;
  G4double zmin = -fScint_z/2. - dz/2.;
  G4int k=0;
 
  z = -fScint_z/2. - height_pmt;      //front
  PlacePMTs(fPmt_log,0,x,y,dx,dy,xmin,ymin,fNx,fNy,x,y,z,k);

  G4RotationMatrix* rm_z = new G4RotationMatrix();
  rm_z->rotateY(180*deg);
  z = fScint_z/2. + height_pmt;       //back
  PlacePMTs(fPmt_log,rm_z,x,y,dx,dy,xmin,ymin,fNx,fNy,x,y,z,k);
 
  G4RotationMatrix* rm_y1 = new G4RotationMatrix();
  rm_y1->rotateY(-90*deg);
  x = -fScint_x/2. - height_pmt;      //left
  PlacePMTs(fPmt_log,rm_y1,y,z,dy,dz,ymin,zmin,fNy,fNz,x,y,z,k);

  G4RotationMatrix* rm_y2 = new G4RotationMatrix();
  rm_y2->rotateY(90*deg);
  x = fScint_x/2. + height_pmt;      //right
  PlacePMTs(fPmt_log,rm_y2,y,z,dy,dz,ymin,zmin,fNy,fNz,x,y,z,k);
 
  G4RotationMatrix* rm_x1 = new G4RotationMatrix();
  rm_x1->rotateX(90*deg);
  y = -fScint_y/2. - height_pmt;     //bottom
  PlacePMTs(fPmt_log,rm_x1,x,z,dx,dz,xmin,zmin,fNx,fNz,x,y,z,k);

  G4RotationMatrix* rm_x2 = new G4RotationMatrix();
  rm_x2->rotateX(-90*deg);
  y = fScint_y/2. + height_pmt;      //top
  PlacePMTs(fPmt_log,rm_x2,x,z,dx,dz,xmin,zmin,fNx,fNz,x,y,z,k);
 
  VisAttributes();
  SurfaceProperties();

  SetLogicalVolume(fHousing_log);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume::CopyValues(){
  fScint_x=fConstructor->GetScintX();
  fScint_y=fConstructor->GetScintY();
  fScint_z=fConstructor->GetScintZ();
  fD_mtl=fConstructor->GetHousingThickness();
  fNx=fConstructor->GetNX();
  fNy=fConstructor->GetNY();
  fNz=fConstructor->GetNZ();
  fOuterRadius_pmt=fConstructor->GetPMTRadius();
  fSphereOn=fConstructor->GetSphereOn();
  fRefl=fConstructor->GetHousingReflectivity();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume::PlacePMTs(G4LogicalVolume* pmt_log,
                              G4RotationMatrix *rot,
                              G4double &a, G4double &b, G4double da,
                              G4double db, G4double amin,
                              G4double bmin, G4int na, G4int nb,
                              G4double &x, G4double &y, G4double &z,
                              G4int &k){
/*PlacePMTs : a different way to parameterize placement that does not depend on
  calculating the position from the copy number

  pmt_log = logical volume for pmts to be placed
  rot = rotation matrix to apply
  a,b = coordinates to vary(ie. if varying in the xy plane then pass x,y)
  da,db = value to increment a,b by
  amin,bmin = start values for a,b
  na,nb = number of repitions in a and b
  x,y,z = just pass x,y, and z by reference (the same ones passed for a,b)
  k = copy number to start with
  sd = sensitive detector for pmts
*/
  a=amin;
  for(G4int j=1;j<=na;j++){
    a+=da;
    b=bmin;
    for(G4int i=1;i<=nb;i++){
      b+=db;
      new G4PVPlacement(rot,G4ThreeVector(x,y,z),pmt_log,"pmt",
                        fHousing_log,false,k);
      fPmtPositions.push_back(G4ThreeVector(x,y,z));
      k++;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume::VisAttributes(){
  G4VisAttributes* housing_va = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
  fHousing_log->SetVisAttributes(housing_va);

  G4VisAttributes* sphere_va = new G4VisAttributes();
  sphere_va->SetForceSolid(true);
  fSphere_log->SetVisAttributes(sphere_va);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume::SurfaceProperties(){
  G4double ephoton[] = {7.0*eV, 7.14*eV};
  const G4int num = sizeof(ephoton)/sizeof(G4double);

  //**Scintillator housing properties
  G4double reflectivity[] = {fRefl, fRefl};
  assert(sizeof(reflectivity) == sizeof(ephoton));
  G4double efficiency[] = {0.0, 0.0};
  assert(sizeof(efficiency) == sizeof(ephoton));
  G4MaterialPropertiesTable* scintHsngPT = new G4MaterialPropertiesTable();
  scintHsngPT->AddProperty("REFLECTIVITY", ephoton, reflectivity, num);
  scintHsngPT->AddProperty("EFFICIENCY", ephoton, efficiency, num);
  G4OpticalSurface* OpScintHousingSurface =
    new G4OpticalSurface("HousingSurface",unified,polished,dielectric_metal);
  OpScintHousingSurface->SetMaterialPropertiesTable(scintHsngPT);
 
  //**Sphere surface properties
  G4double sphereReflectivity[] = {1.0, 1.0};
  assert(sizeof(sphereReflectivity) == sizeof(ephoton));
  G4double sphereEfficiency[] = {0.0, 0.0};
  assert(sizeof(sphereEfficiency) == sizeof(ephoton));
  G4MaterialPropertiesTable* spherePT = new G4MaterialPropertiesTable();
  spherePT->AddProperty("REFLECTIVITY", ephoton, sphereReflectivity, num);
  spherePT->AddProperty("EFFICIENCY", ephoton, sphereEfficiency, num);
  G4OpticalSurface* OpSphereSurface =
    new G4OpticalSurface("SphereSurface",unified,polished,dielectric_metal);
  OpSphereSurface->SetMaterialPropertiesTable(spherePT);
 
  //**Photocathode surface properties
  G4double photocath_EFF[]={1.,1.}; //Enables 'detection' of photons
  assert(sizeof(photocath_EFF) == sizeof(ephoton));
  G4double photocath_ReR[]={1.92,1.92};
  assert(sizeof(photocath_ReR) == sizeof(ephoton));
  G4double photocath_ImR[]={1.69,1.69};
  assert(sizeof(photocath_ImR) == sizeof(ephoton));
  G4MaterialPropertiesTable* photocath_mt = new G4MaterialPropertiesTable();
  photocath_mt->AddProperty("EFFICIENCY",ephoton,photocath_EFF,num);
  photocath_mt->AddProperty("REALRINDEX",ephoton,photocath_ReR,num);
  photocath_mt->AddProperty("IMAGINARYRINDEX",ephoton,photocath_ImR,num);
  G4OpticalSurface* photocath_opsurf=
    new G4OpticalSurface("photocath_opsurf",glisur,polished,
                         dielectric_metal);
  photocath_opsurf->SetMaterialPropertiesTable(photocath_mt);

  //**Create logical skin surfaces
  new G4LogicalSkinSurface("photocath_surf",fHousing_log,
                           OpScintHousingSurface);
  new G4LogicalSkinSurface("sphere_surface",fSphere_log,OpSphereSurface);
  new G4LogicalSkinSurface("photocath_surf",fPhotocath_log,photocath_opsurf);
}
