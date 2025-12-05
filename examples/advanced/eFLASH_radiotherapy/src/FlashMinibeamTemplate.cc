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
// This is the first version of Flash, a Geant4-based application
//
//
//////////////////////////////////////////////////////////////////////////////////////////////

#include "FlashMinibeamTemplate.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4NistElementBuilder.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4RunManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4MultiUnion.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"
#include "G4UserLimits.hh"
#include "G4Region.hh"
#include <cmath>
FlashMinibeamTemplate::FlashMinibeamTemplate(G4VPhysicalVolume *physicalTreatmentRoom, G4double x, G4double r)
  : hight(0.),motherPhys(physicalTreatmentRoom),X_center(x), y_center(0.),z_center(0.),d_between_holes(0.),hole_side(0.),radius(r),field_side(0)
       {
  ConstructCollimator();
}

FlashMinibeamTemplate::~FlashMinibeamTemplate() {}

void FlashMinibeamTemplate::ConstructCollimator() {
  // Sets default geometry and materials
  SetMaterial();
  //ConstructColl_full(); //FULL
  // ConstructColl_template(); //GRID
  // Construct_hole();  //SINGLE HOLE
  ConstructColl_template_planar(); //PLANAR
}

void FlashMinibeamTemplate::SetMaterial() {
  //X_center = 0 *cm; 
  y_center = 0.*cm;
  z_center = 0.*cm;
  d_between_holes=2; //mm
  hole_side = 1; //mm
  maxStep = 0.1 * mm;
  field_side=10; //mm
  TemplateRegion = new G4Region("Template_region");
  fStepLimit = new G4UserLimits(maxStep);
  blue = new G4VisAttributes(G4Colour(0., 0., 1.));
  red = new G4VisAttributes(G4Colour(0., 1., 0.));
  blue->SetVisibility(true);
  red->SetVisibility(true);
  G4NistManager *man = G4NistManager::Instance();
  G4bool isotopes = false;
  TEFLON = man->FindOrBuildMaterial("G4_TEFLON", isotopes);
  TUNGSTEN = man->FindOrBuildMaterial("G4_W", isotopes);
}

void FlashMinibeamTemplate::ConstructColl_full() {
G4double phi1 = 90. * deg;
  G4RotationMatrix rm1;
  rm1.rotateY(phi1);
  const G4double outerRad = radius;
  const G4double innRadius = 0. * mm;
	hight = 5 * cm;
  const G4double startAngle= 0. * deg;
  const G4double spanningAngle = 360. * deg;
  const G4double XPosition = X_center+hight;
  G4VSolid * solid =
      new G4Tubs("Solid_template", innRadius, outerRad, hight,
                 startAngle, spanningAngle);
  G4LogicalVolume *log = new G4LogicalVolume(
      solid, TUNGSTEN, "Log_template", 0, 0, 0);
   new G4PVPlacement(
      G4Transform3D(rm1, G4ThreeVector((XPosition), 0., 0.)), "physi_template",
      log, motherPhys, false, 0);
  log->SetVisAttributes(blue);
  log->SetRegion(TemplateRegion);
  TemplateRegion->AddRootLogicalVolume(log);
  log->SetUserLimits(fStepLimit);
}

void FlashMinibeamTemplate::ConstructColl_template() {  //GRID
G4double phi2 = 90. * deg;
  G4RotationMatrix rm2;
  rm2.rotateY(phi2);
	hight = 0.25 * cm;
  const G4double XPosition = X_center+hight;
G4Box *t1 = new G4Box("t1_", radius, radius, hight); 
G4Box *Air_box = new G4Box("air_templ",hole_side/2,hole_side/2,hight+0.1*mm);
G4MultiUnion* munion_solid = new G4MultiUnion("Boxes_Union");
G4int n_holes = (G4int) round(field_side/(d_between_holes));
G4RotationMatrix rotm = G4RotationMatrix();
G4cout<< "NUMBER:"<<n_holes<<G4endl;
for (G4int i= 0; i<n_holes;i++){
  for (G4int j=0;j<n_holes;j++){
    G4double xi = (-(n_holes -1)/2+ i)*(d_between_holes) ;
    G4double yi = ((n_holes -1)/2- j)*(d_between_holes) ;
G4ThreeVector position = G4ThreeVector(xi*mm,yi*mm,0.);
G4Transform3D tr = G4Transform3D(rotm,position);
munion_solid->AddNode(*Air_box,tr);
}
 }
munion_solid->Voxelize();
// Associate it to a logical volume as a normal solid
//
G4RotationMatrix rotm_t2 = G4RotationMatrix();
    rotm_t2.rotateX(0 * deg);
    G4ThreeVector zTrans(0, 0, 0);
    G4SubtractionSolid *hollowcover =
        new G4SubtractionSolid("template", t1, munion_solid, 0, zTrans);
G4LogicalVolume *logic =
        new G4LogicalVolume(hollowcover, TUNGSTEN, "Log_template", 0, 0, 0);
 new G4PVPlacement(
      G4Transform3D(rm2, G4ThreeVector((XPosition), 0., 0.)), "cover1phys",
      logic, motherPhys, false, 0);
      logic->SetVisAttributes(blue);
      logic->SetRegion(TemplateRegion);
      TemplateRegion->AddRootLogicalVolume(logic);
      logic->SetUserLimits(fStepLimit);
}
void FlashMinibeamTemplate::Construct_hole(){ //Single hole
G4double phi2 = 90. * deg;
  G4RotationMatrix rm2;
  rm2.rotateY(phi2);
  const G4double outerRad = radius;
  const G4double innRadius = 0. * mm;
	hight = 0.25 * cm;
  const G4double startAngle= 0. * deg;
  const G4double spanningAngle = 360. * deg;
  const G4double XPosition = X_center+hight;    
  G4VSolid *t1 = new G4Tubs("t1_", innRadius, outerRad, hight,
                 startAngle, spanningAngle);
G4Box *Air_box = new G4Box("air_templ",0.5*cm,0.5*cm,hight+1*mm);
G4RotationMatrix rotm_t2 = G4RotationMatrix();
    rotm_t2.rotateX(0 * deg);
    //G4ThreeVector zTrans(lato/2, -lato/2, 0);
    G4ThreeVector zTrans(0, 0, 0);
    G4SubtractionSolid *hollowcover =
        new G4SubtractionSolid("template", t1, Air_box, 0, zTrans);
G4LogicalVolume *logic =
        new G4LogicalVolume(hollowcover, TUNGSTEN, "Log_template", 0, 0, 0);
    new G4PVPlacement(
      G4Transform3D(rm2, G4ThreeVector((XPosition), 0., 0.)), "cover1phys",
      logic, motherPhys, false, 0);
        logic->SetVisAttributes(blue);
  logic->SetRegion(TemplateRegion);
  TemplateRegion->AddRootLogicalVolume(logic);
  logic->SetUserLimits(fStepLimit);
}

void FlashMinibeamTemplate::ConstructColl_template_planar() {
G4double phi2 = 90. * deg;
  G4RotationMatrix rm2;
  rm2.rotateY(phi2);
	hight = 0.25 * cm;
  const G4double XPosition = X_center+hight;
G4Box *t1 = new G4Box("t1_", radius, radius, hight);
 G4Box *Air_box = new G4Box("air_templ",hole_side/2,field_side/2*mm,hight+0.1*mm);
G4MultiUnion* munion_solid = new G4MultiUnion("Boxes_Union");
G4int n_holes = (G4int) round(field_side/(d_between_holes))/2; //TOTAL HOLE NUMBERS=n_holes*2+1
G4RotationMatrix rotm = G4RotationMatrix();
G4ThreeVector position = G4ThreeVector(0.,0.,0.);
G4Transform3D tr = G4Transform3D(rotm,position);
munion_solid->AddNode(*Air_box,tr);
  for (G4int j=1;j<n_holes;j++){
G4double xi = j*(d_between_holes) ;
//G4double yi = ((n_holes -1)/2- j)*(d_between_holes) ;
 position = G4ThreeVector(xi*mm,0.,0.);
 tr = G4Transform3D(rotm,position);
munion_solid->AddNode(*Air_box,tr);
 position = G4ThreeVector(-xi*mm,0.,0.);
 tr = G4Transform3D(rotm,position);
munion_solid->AddNode(*Air_box,tr);
}
munion_solid->Voxelize();
// Associate it to a logical volume as a normal solid
//
G4RotationMatrix rotm_t2 = G4RotationMatrix();
    rotm_t2.rotateX(0 * deg);
    //G4ThreeVector zTrans(lato/2, -lato/2, 0);
    G4ThreeVector zTrans(0, 0, 0);
    G4SubtractionSolid *hollowcover =
        new G4SubtractionSolid("template", t1, munion_solid, 0, zTrans);
G4LogicalVolume *logic =
        new G4LogicalVolume(hollowcover, TUNGSTEN, "Log_template", 0, 0, 0);
    new G4PVPlacement(
      G4Transform3D(rm2, G4ThreeVector((XPosition), 0., 0.)), "cover1phys",
      logic, motherPhys, false, 0);
        logic->SetVisAttributes(blue);
  logic->SetRegion(TemplateRegion);
  TemplateRegion->AddRootLogicalVolume(logic);
  logic->SetUserLimits(fStepLimit);
}

