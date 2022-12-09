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

#include "Applicator.hh"
#include "FlashDetectorConstruction.hh"
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
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"

Applicator::Applicator(G4VPhysicalVolume *physicalTreatmentRoom)
    : fMotherPhys(physicalTreatmentRoom),

      
      solidFTFlash(0), physiFTFlash(0) {
  ConstructCollimator(fMotherPhys);


  
}

Applicator::~Applicator() {}

void Applicator::ConstructCollimator(G4VPhysicalVolume *) {
  // Sets default geometry and materials
  SetDefaultDimensions();
   SetOuterRadius(55*mm);
   SetApplicatorLength(300*mm);
  // Construct the whole Applicator Beam Line
  ConstructApplicator();
}

void Applicator::SetDefaultDimensions() {

  white = new G4VisAttributes(G4Colour());
  white->SetVisibility(true);

  blue = new G4VisAttributes(G4Colour(0., 0., 1.));
  blue->SetVisibility(true);

  gray = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
  gray->SetVisibility(true);

  red = new G4VisAttributes(G4Colour(1., 0., 0.));
  red->SetVisibility(true);

  yellow = new G4VisAttributes(G4Colour(1., 1., 0.));
  yellow->SetVisibility(true);

  green = new G4VisAttributes(G4Colour(25 / 255., 255 / 255., 25 / 255.));
  green->SetVisibility(true);

  darkGreen = new G4VisAttributes(G4Colour(0 / 255., 100 / 255., 0 / 255.));
  darkGreen->SetVisibility(true);

  darkOrange3 =
      new G4VisAttributes(G4Colour(205 / 255., 102 / 255., 000 / 255.));
  darkOrange3->SetVisibility(true);

  skyBlue = new G4VisAttributes(G4Colour(135 / 255., 206 / 255., 235 / 255.));
  skyBlue->SetVisibility(true);

  magenta = new G4VisAttributes(G4Colour(255 / 255., 0 / 255., 255 / 255.));
  magenta->SetVisibility(true);

  
  fInitial_pos = -100*cm; //set the same position in FlashPrimaryGeneratorAction.cc
  // Geometry  APPLICATOR DEFAULTS

 

  G4double defaultinnerRadiusFirstApplicatorFlash =
      fOuterRadiusFirstApplicatorFlash - 5. * mm;
   fInnerRadiusFirstApplicatorFlash = defaultinnerRadiusFirstApplicatorFlash;

  // DEFAULT DEFINITION OF THE MATERIALS

  // ELEMENTS
  G4double density;
  G4int ncomponents;
  G4bool isotopes = false;
  aluminumNist =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_Al", isotopes);


  Fe = G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe");
    
  PVDF = new G4Material("PVDF",  density=1780 *kg/m3,  ncomponents=3);
  
  
  
  PVDF->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 34 * perCent);
  PVDF->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), 33 * perCent);
  PVDF->AddElement(G4NistManager::Instance()->FindOrBuildElement("F"), 33 * perCent);

  
  
  FILM= new G4Material("FILM", density=1430 *kg/m3, ncomponents = 4);
	
  FILM->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 69 * perCent);
  FILM->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), 3 * perCent);
  FILM->AddElement(G4NistManager::Instance()->FindOrBuildElement("N"), 7 * perCent);
  FILM->AddElement(G4NistManager::Instance()->FindOrBuildElement("O"), 21 * perCent);
	



  G4Material *galacticNist =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic", isotopes);
  PMMA =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS", isotopes);

  G4Material *titanioNist =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_Ti", isotopes);


  // MATERIAL ASSIGNMENT

  // Material of the APPLICATOR Flash
  fFirstApplicatorMaterialFlash = PMMA;
 
  // Titanium window
  FTFlashMaterialFlash = titanioNist;

  // Vacuum Source
  VSFlashMaterialFlash = galacticNist;
}

void Applicator::ConstructApplicator() {

  // Components of the Applicator

  FlashBeamLineVacuumSource();
  FlashBeamLineTitaniumWindows();
  FlashVWAlcover();
  FlashAlCover2();
  FlashExitBit();
  FlashToroid();
  OverCover();

  MonitorChamber();
  Flash_connector();
  Bigconnector();
  Bigconnector2();
  FlashBeamLineApplicator();//modify this function for applicator lenght


 
}

void Applicator::FlashBeamLineVacuumSource() {
  // ---------------------------------------------------------------//
  //                     Vacuum Source                             //
  // ---------------------------------------------------------------//

  G4double phi1 = 90. * deg;

  G4RotationMatrix rm1;
  rm1.rotateY(phi1);

  fOutRadiusVSFlash = 20 * mm;
  const G4double innRadiusVSFlash = 0. * mm;
   fHightVSFlash = 8 * mm;
  const G4double startAngleVSFlash = 0. * deg;
  const G4double spanningAngleVSFlash = 360. * deg;
  fXPositionVSFlash = fInitial_pos-fHightVSFlash-0.055/2*mm; //0.055 is titanium window

  solidVSFlash =
      new G4Tubs("VSFlash", innRadiusVSFlash, fOutRadiusVSFlash, fHightVSFlash,
                 startAngleVSFlash, spanningAngleVSFlash);

  G4LogicalVolume *logVSFlash = new G4LogicalVolume(
      solidVSFlash, VSFlashMaterialFlash, "VSFlash", 0, 0, 0);

  physiVSFlash = new G4PVPlacement(
      G4Transform3D(rm1, G4ThreeVector((fXPositionVSFlash), 0., 0.)), "VSFlash",
      logVSFlash, fMotherPhys, false, 0);

  logVSFlash->SetVisAttributes(green);
}

void Applicator::FlashBeamLineTitaniumWindows() {
  // ---------------------------------------------------------------//
  //                     Titanium Window                        //
  // ---------------------------------------------------------------//
  // with just this piece ssd=1.6cm
  G4double phi2 = 90. * deg;

  G4RotationMatrix rm2;
  rm2.rotateY(phi2);

   fOutRadiusFTFlash = fOutRadiusVSFlash;
  const G4double innRadiusFTFlash = 19 * mm;
   fHightFTFlash = 0.055/2 * mm;
  const G4double startAngleFTFlash = 0. * deg;
  const G4double spanningAngleFTFlash = 360. * deg;
  const G4double XPositionFTFlash = fInitial_pos ;

  solidFTFlash =
      new G4Tubs("FTFlash", innRadiusFTFlash, fOutRadiusFTFlash, fHightFTFlash,
                 startAngleFTFlash, spanningAngleFTFlash);

  G4LogicalVolume *logFTFlash = new G4LogicalVolume(
      solidFTFlash, FTFlashMaterialFlash, "FTFlash", 0, 0, 0);

  physiFTFlash = new G4PVPlacement(
      G4Transform3D(rm2, G4ThreeVector((XPositionFTFlash), 0., 0.)), "FTFlash",
      logFTFlash, fMotherPhys, false, 0);

  logFTFlash->SetVisAttributes(yellow);
}
void Applicator::FlashVWAlcover(){

 G4double phi2 = 90. * deg;

  G4RotationMatrix rm2;
  rm2.rotateY(phi2);


  const G4double innRadius = fOutRadiusVSFlash;
    fOutRadius = innRadius+8*mm;
  const G4double hight = fHightVSFlash;
  const G4double startAngle = 0. * deg;
  const G4double spanningAngle = 360. * deg;
  const G4double XPosition = fXPositionVSFlash ;

   G4VSolid * solid =
      new G4Tubs("cover1", innRadius, fOutRadius, hight,
                 startAngle, spanningAngle);

  G4LogicalVolume *log = new G4LogicalVolume(
      solid, aluminumNist, "cover1log", 0, 0, 0);

 new G4PVPlacement(
      G4Transform3D(rm2, G4ThreeVector((XPosition), 0., 0.)), "cover1phys",
      log, fMotherPhys, false, 0);

  log->SetVisAttributes(white);


}
void Applicator::FlashAlCover2(){

G4double phi2 = 90. * deg;

  G4RotationMatrix rm2;
  rm2.rotateY(phi2);

  
  const G4double innRadius = fOutRadiusVSFlash;

  const G4double hight = fHightVSFlash+fHightFTFlash;
  const G4double startAngle = 0. * deg;
  const G4double spanningAngle = 360. * deg;
  const G4double XPosition = fInitial_pos+hight+fHightFTFlash;

   G4VSolid * solid =
      new G4Tubs("cover1", innRadius, fOutRadius, hight,
                 startAngle, spanningAngle);

  G4LogicalVolume *log = new G4LogicalVolume(
      solid, aluminumNist, "cover1log", 0, 0, 0);

new G4PVPlacement(
      G4Transform3D(rm2, G4ThreeVector((XPosition), 0., 0.)), "cover1phys",
      log, fMotherPhys, false, 0);

  log->SetVisAttributes(red);
  
  fInitial_pos=fInitial_pos + fHightFTFlash;
}
void Applicator::FlashExitBit(){

G4double phi2 = 90. * deg;

  G4RotationMatrix rm2;
  rm2.rotateY(phi2);

  
  const G4double innRadius = 0*mm;
  fOutRadius = fOutRadiusVSFlash;
  const G4double hight = 16/2*mm;
  const G4double startAngle = 0. * deg;
  const G4double spanningAngle = 360. * deg;
  const G4double XPosition = fInitial_pos+hight;


    G4VSolid *t1 = new G4Tubs("t1", innRadius, fOutRadius, hight,
                 startAngle, spanningAngle);
                             
    G4VSolid *t2 = new G4Cons("t2", 0*mm,13/2*mm, 0*mm,38/2*mm,16.1/2*mm, startAngle,spanningAngle);
    
    G4RotationMatrix rotm_t2 = G4RotationMatrix();
    rotm_t2.rotateX(0 * deg);
    G4ThreeVector zTrans(0, 0, 0);

    G4SubtractionSolid *hollowcover =
        new G4SubtractionSolid("hollowcover_log", t1, t2, 0, zTrans);

    
    G4LogicalVolume *logic =
        new G4LogicalVolume(hollowcover, aluminumNist, "hollowcover", 0, 0, 0);
        
    new G4PVPlacement(
      G4Transform3D(rm2, G4ThreeVector((XPosition), 0., 0.)), "cover1phys",
      logic, fMotherPhys, false, 0);

  logic->SetVisAttributes(darkOrange3);
fInitial_pos=XPosition+hight;
}
void Applicator::FlashToroid(){

G4double phi2 = 90. * deg;

  G4RotationMatrix rm2;
  rm2.rotateY(phi2);

  
  const G4double innRadius = 50.8/2*mm;
 fToroid_outRadius = innRadius + 45*mm;
  fToroid_hight = 50.8/2*mm;
  const G4double startAngle = 0. * deg;
  const G4double spanningAngle = 360. * deg;
  fToroid_XPosition = fInitial_pos+fToroid_hight + 8*mm;

   G4VSolid * solid =
      new G4Tubs("toroid", innRadius, fToroid_outRadius, fToroid_hight,
                 startAngle, spanningAngle);

  G4LogicalVolume *log = new G4LogicalVolume(
      solid, Fe, "toroidlog", 0, 0, 0);

  new G4PVPlacement(
      G4Transform3D(rm2, G4ThreeVector((fToroid_XPosition), 0., 0.)), "toroidphys",
      log, fMotherPhys, false, 0);

  log->SetVisAttributes(blue);

fInitial_pos=fToroid_XPosition+fToroid_hight;
}
void Applicator::OverCover(){

G4double phi2 = 90. * deg;

  G4RotationMatrix rm2;
  rm2.rotateY(phi2);

  
  const G4double innRadius = fOutRadius;
  const G4double out_Radius = 5*fOutRadius ;
  fBigcover_hight = 7.5*mm;
  const G4double startAngle = 0. * deg;
  const G4double spanningAngle = 360. * deg;
   fBigcover_XPosition = fXPositionVSFlash - fHightVSFlash+fBigcover_hight;

   G4VSolid * solid =
      new G4Tubs("coverbig", innRadius, out_Radius, fBigcover_hight,
                 startAngle, spanningAngle);

  G4LogicalVolume *log = new G4LogicalVolume(
      solid, aluminumNist, "coverbig_log", 0, 0, 0);

  new G4PVPlacement(
      G4Transform3D(rm2, G4ThreeVector((fBigcover_XPosition), 0., 0.)), "coverbig_phys",
      log, fMotherPhys, false, 0);

  log->SetVisAttributes(skyBlue);
  
  
  const G4double innRadius_2 = fToroid_outRadius;
  const G4double out_Radius_2 = innRadius_2 + 1.2*cm ;
  const G4double fBigcover_hight_2 = 30*mm;
  const G4double startAngle_2 = 0. * deg;
  const G4double spanningAngle_2 = 360. * deg;
   const double fBigcover_XPosition_2 = fBigcover_XPosition+fBigcover_hight_2+fBigcover_hight;

   G4VSolid * solid_2 =
      new G4Tubs("coverbig_2", innRadius_2, out_Radius_2, fBigcover_hight_2,
                 startAngle_2, spanningAngle_2);

  G4LogicalVolume *log_2 = new G4LogicalVolume(
      solid_2, aluminumNist, "coverbig_log", 0, 0, 0);

  new G4PVPlacement(
      G4Transform3D(rm2, G4ThreeVector((fBigcover_XPosition_2), 0., 0.)), "coverbig_phys",
      log_2, fMotherPhys, false, 0);

  log_2->SetVisAttributes(green);
  
  
}

void Applicator::OverCover2() {

G4double phi2 = 90. * deg;

  G4RotationMatrix rm2;
  rm2.rotateY(phi2);

  
  const G4double innRadius = fToroid_outRadius;
  const G4double out_Radius = innRadius+40*mm ;
  const G4double hight = 34*mm;
  const G4double startAngle = 0. * deg;
  const G4double spanningAngle = 360. * deg;
  const G4double XPosition = fBigcover_XPosition+fBigcover_hight+hight;

   G4VSolid * solid =
      new G4Tubs("coverbig", innRadius, out_Radius, hight,
                 startAngle, spanningAngle);

  G4LogicalVolume *log = new G4LogicalVolume(
      solid, aluminumNist, "coverbig_log", 0, 0, 0);

new G4PVPlacement(
      G4Transform3D(rm2, G4ThreeVector((XPosition), 0., 0.)), "coverbig_phys",
      log, fMotherPhys, false, 0);

  log->SetVisAttributes(yellow);

}

void Applicator::MonitorChamber(){

G4double phi2 = 90. * deg;

  G4RotationMatrix rm2;
  rm2.rotateY(phi2);

  
  const G4double innRadius = 20*mm;
  const G4double out_Radius = innRadius+1.7*mm ;
  const G4double hight = 3*mm;
  const G4double startAngle = 0. * deg;
  const G4double spanningAngle = 360. * deg;
   G4double XPosition = fInitial_pos+hight;

   G4VSolid * solid =
      new G4Tubs("first", innRadius, out_Radius, hight,
                 startAngle, spanningAngle);



G4LogicalVolume *log = new G4LogicalVolume(
      solid, PVDF, "chamberfirst_log", 0, 0, 0);

  new G4PVPlacement(
      G4Transform3D(rm2, G4ThreeVector((XPosition), 0., 0.)), "coverbig_phys",
      log, fMotherPhys, false, 0);

  log->SetVisAttributes(red);
  
  G4VSolid * solid_pvdf=
      new G4Tubs("s_pvdf", innRadius, out_Radius, 0.5*mm,
                 startAngle, spanningAngle);
    
      
      G4LogicalVolume *log_pvdf= new G4LogicalVolume(
      solid_pvdf, PVDF, "pvdf_log", 0, 0, 0);
                 
  G4VSolid * solid_film =
      new G4Tubs("s_film", innRadius, out_Radius, 0.05/2*mm,
                 startAngle, spanningAngle);
   G4LogicalVolume *log_film = new G4LogicalVolume(
      solid_film, FILM, "ka_log", 0, 0, 0);
      
  G4VSolid * solid_al =
      new G4Tubs("s_al", innRadius, out_Radius, 0.005/2*mm,
                 startAngle, spanningAngle);          
  
  G4LogicalVolume *log_al = new G4LogicalVolume(
      solid_al, aluminumNist, "al_log", 0, 0, 0);
  XPosition=XPosition+hight;

  G4int j=0;
  for(G4int i = 0;i<3;i++){
  XPosition=XPosition+(i+1)*0.05/2*mm;
  new G4PVPlacement(
      G4Transform3D(rm2, G4ThreeVector((XPosition), 0., 0.)), "ka_phys",
      log_film, fMotherPhys, false, j);
  XPosition=XPosition+(i+1)*0.005/2*mm;
  new G4PVPlacement(
      G4Transform3D(rm2, G4ThreeVector((XPosition), 0., 0.)), "al_phys",
      log_al, fMotherPhys, false, j);
      
  XPosition=XPosition+(i+1)*1/2*mm;
  new G4PVPlacement(
      G4Transform3D(rm2, G4ThreeVector((XPosition), 0., 0.)), "pvdf_phys",
      log_pvdf, fMotherPhys, false, j);
  
  
  
  }
  fChamberpos = XPosition +1/2*mm;
  log_film->SetVisAttributes(green);
    log_al->SetVisAttributes(blue);
      log_pvdf->SetVisAttributes(yellow);


}

void Applicator::Flash_connector(){

G4double phi2 = 90. * deg;

  G4RotationMatrix rm2;
  rm2.rotateY(phi2);

  
  const G4double innRadius = 10*cm;
  const G4double out_Radius = innRadius+2.5*cm ;
  const G4double hight = 15*mm;
  const G4double startAngle = 0. * deg;
  const G4double spanningAngle = 360. * deg;
  const G4double XPosition = fChamberpos+hight;
  
  G4VSolid * solid =
      new G4Tubs("cover", innRadius, out_Radius, hight,
                 startAngle, spanningAngle);

  G4LogicalVolume *log = new G4LogicalVolume(
      solid, aluminumNist, "coverbig_log", 0, 0, 0);

 new G4PVPlacement(
      G4Transform3D(rm2, G4ThreeVector((XPosition), 0., 0.)), "coverbig_phys",
      log, fMotherPhys, false, 0);

  log->SetVisAttributes(magenta);
  
  
  G4VSolid * solid_ =
      new G4Tubs("littlecover", 22*mm, innRadius, 0.2*mm,
                 startAngle, spanningAngle);
  G4LogicalVolume *log_ = new G4LogicalVolume(
      solid_, aluminumNist, "covers_log", 0, 0, 0);
      
    new G4PVPlacement(
      G4Transform3D(rm2, G4ThreeVector((fChamberpos+0.2*mm), 0., 0.)), "coverl_phys",
      log_, fMotherPhys, false, 0);  
log_->SetVisAttributes(green);

fInitial_pos=XPosition+hight;
}

void Applicator::Bigconnector() {

G4double phi2 = 90. * deg;

  G4RotationMatrix rm2;
  rm2.rotateY(phi2);

  
  const G4double innRadius = 10*cm;
  const G4double out_Radius = innRadius+30*mm ;
  const G4double hight = 7.05*cm + 0.0075/2*mm;
  const G4double startAngle = 0. * deg;
  const G4double spanningAngle = 360. * deg;
  const G4double XPosition = fInitial_pos+hight;

   G4VSolid * solid =
      new G4Tubs("coverbig", innRadius, out_Radius, hight,
                 startAngle, spanningAngle);

  G4LogicalVolume *log = new G4LogicalVolume(
      solid, aluminumNist, "coverbig_log", 0, 0, 0);

 new G4PVPlacement(
      G4Transform3D(rm2, G4ThreeVector((XPosition), 0., 0.)), "coverbig_phys",
      log, fMotherPhys, false, 0);

  log->SetVisAttributes(red);
fInitial_pos=XPosition+hight;
}
void Applicator::Bigconnector2() {

G4double phi2 = 90. * deg;

  G4RotationMatrix rm2;
  rm2.rotateY(phi2);

  
  const G4double innRadius = 6*cm;
  const G4double out_Radius = innRadius+70*mm ;
  const G4double hight = 4.4*cm-12/4*mm;
  const G4double startAngle = 0. * deg;
  const G4double spanningAngle = 360. * deg;
  const G4double XPosition = fInitial_pos+hight;

   G4VSolid * solid =
      new G4Tubs("coverbig", innRadius, out_Radius, hight,
                 startAngle, spanningAngle);



      
      G4LogicalVolume *log = new G4LogicalVolume(
      solid, PVDF, "coverbig_log", 0, 0, 0);

 new G4PVPlacement(
      G4Transform3D(rm2, G4ThreeVector((XPosition), 0., 0.)), "coverbig_phys",
      log, fMotherPhys, false, 0);

  log->SetVisAttributes(blue);
fInitial_pos=XPosition+hight;
}

void Applicator::Bigconnector3() {//this is a piece of the applicator that connects the tube to the optics

G4double phi2 = 90. * deg;

  G4RotationMatrix rm2;
  rm2.rotateY(phi2);

  
  const G4double innRadius = 6*cm;
  const G4double out_Radius = innRadius+60*mm ;
  const G4double hight = 3.4*cm-11/4*mm;
  const G4double startAngle = 0. * deg;
  const G4double spanningAngle = 360. * deg;
  const G4double XPosition = fInitial_pos+hight;


                 
  G4VSolid *t1 = new G4Tubs("t1_", 0*mm, out_Radius, hight,
                 startAngle, spanningAngle);
                             
    G4VSolid *t2 = new G4Cons("t2_", 0*mm,60*mm, 0*mm,fInnerRadiusFirstApplicatorFlash,hight +0.1*mm, startAngle,spanningAngle);
    
   
    G4ThreeVector zTrans(0, 0, 0);

    G4SubtractionSolid *hollowcover =
        new G4SubtractionSolid("hollowcover_log_", t1, t2, 0, zTrans);

  G4LogicalVolume *log = new G4LogicalVolume(
      hollowcover, PMMA, "coverbig_log_", 0, 0, 0);

 new G4PVPlacement(
      G4Transform3D(rm2, G4ThreeVector((XPosition), 0., 0.)), "coverbig_phys_",
      log, fMotherPhys, false, 0);

  log->SetVisAttributes(yellow);
fInitial_pos=XPosition+hight;
}

void Applicator::FlashBeamLineApplicator() {

  


  // hightFinalApplicatorFlash = 300 * mm;//modify this for length of applicator
   if (fHightFinalApplicatorFlash != 0*mm){//set to zero if you do not want the applicator
     Bigconnector3(); //comment this line to remove the applicator piece
  const G4double startAngleFirstApplicatorFlash = 0. * deg;
  const G4double spanningAngleFirstApplicatorFlash = 360. * deg;
   fFinalApplicatorXPositionFlash = fInitial_pos+fHightFinalApplicatorFlash;

  G4double phi6 = 90. * deg;

  G4RotationMatrix rm6;
  rm6.rotateY(phi6);

  fSolidFirstApplicatorFlash = new G4Tubs(
      "FirstApplicatorFlash", fInnerRadiusFirstApplicatorFlash,
      fOuterRadiusFirstApplicatorFlash, fHightFinalApplicatorFlash,
      startAngleFirstApplicatorFlash, spanningAngleFirstApplicatorFlash);

  G4LogicalVolume *logFirstApplicatorFlash = new G4LogicalVolume(
      fSolidFirstApplicatorFlash, fFirstApplicatorMaterialFlash,
      "FirstApplicatorFlash", 0, 0, 0);

  fPhysiFirstApplicatorFlash = new G4PVPlacement(
      G4Transform3D(rm6,
                    G4ThreeVector((fFinalApplicatorXPositionFlash), 0., 0.)),
      "FirstApplicatorFlash", logFirstApplicatorFlash, fMotherPhys, false, 0);

  logFirstApplicatorFlash->SetVisAttributes(magenta); } else{fFinalApplicatorXPositionFlash = fInitial_pos+3*cm;}
} 


void Applicator::SetOuterRadius(G4double radius)
{
 
  fOuterRadiusFirstApplicatorFlash=radius;
  fInnerRadiusFirstApplicatorFlash= fOuterRadiusFirstApplicatorFlash-5*mm;
  
}

void Applicator::SetApplicatorLength(G4double length)
{
  fHightFinalApplicatorFlash=length;
  
}





