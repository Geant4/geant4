//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// --------------------------------------------------------------
//                 GEANT 4 - ULTRA experiment example
// --------------------------------------------------------------
//
// Code developed by:
// B. Tome, M.C. Espirito-Santo, A. Trindade, P. Rodrigues 
//
//    ****************************************************
//    *      UltraDetectorConstruction.cc
//    ****************************************************
//
//    Class used in the definition of the Ultra setup consisting of:
//      - the UVscope detector
//      - an optional reflecting surface
//    Optical photons can reach the UVscope either directly or after reflection in the
//    surface, which can be polished or diffusing.
//    The main part of the UVscope definition is the Fresnel lens construction based
//    on the UltraFresnelLens class.
//
#include <cmath>
#include "UltraDetectorConstruction.hh"
#include "UltraPMTSD.hh"
#include "UltraFresnelLens.hh"

#include "G4SDManager.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UltraDetectorConstruction::UltraDetectorConstruction()
{

 PMTSD   = 0;

 // Sensitive Detector Manager
 SDmanager = G4SDManager::GetSDMpointer();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UltraDetectorConstruction::~UltraDetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* UltraDetectorConstruction::Construct()
{
  ConstructTableMaterials();



//	The experimental Hall
//	---------------------

  G4double World_x = 1.*m;
  G4double World_y = 1.*m;
  G4double World_z = 2*m;

  G4Box * World_box = new G4Box("World",World_x,World_y,World_z);

  // Get Air pointer from static funcion - (G4Material::GetMaterial)

G4String name;
G4Material *Air = G4Material::GetMaterial(name = "Air");
G4LogicalVolume *World_log ;
World_log  = new G4LogicalVolume(World_box,Air,"World",0,0,0);

G4VPhysicalVolume *World_phys ;
World_phys   = new G4PVPlacement(0,G4ThreeVector(),"World",World_log,0,false,0);

   G4VisAttributes* UniverseVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
   UniverseVisAtt->SetVisibility(true);
   UniverseVisAtt->SetForceWireframe(true);
   World_log->SetVisAttributes(UniverseVisAtt);
   World_log->SetVisAttributes (G4VisAttributes::Invisible);



  G4cout << "\n \n \n \n \n \n \n \n \n \n \n \n \n " << G4endl ;

  G4cout << "######################################################" << G4endl ;
  G4cout << "#                                                    #" << G4endl ;
  G4cout << "#                                                    #" << G4endl ;
  G4cout << "#          UltraDetectorConstruction:                #" << G4endl ;
  G4cout << "#                                                    #" << G4endl ;  
  G4cout << "#                                                    #" << G4endl ;  

  G4VPhysicalVolume* chosenVolume;
  chosenVolume = ConstructUVscope(World_phys);


  G4cout << "#                                                    #" << G4endl ;
  G4cout << "#                                                    #" << G4endl ;
  G4cout << "######################################################" << G4endl ;


#ifdef ULTRA_MIRROR_USE

  G4cout << "Using mirror reflecting surface " << G4endl ;

  G4VPhysicalVolume* Mirror ;
  Mirror = ConstructMirror(World_phys);

#elif ULTRA_GROUND_USE

  G4cout << "Using ground reflecting surface " << G4endl ;

  G4VPhysicalVolume* Ground ;
  Ground = ConstructGround(World_phys);

#else

  G4cout << "No reflecting surface used" << G4endl ;

#endif

  return World_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UltraDetectorConstruction::ConstructTableMaterials()
{
  G4double a, z, density;
  G4String name, symbol;
  G4int nel;


//	------------- Elements -------------
  a = 1.01*g/mole;
  G4Element* elH  = new G4Element(name="Hydrogen", symbol="H", z=1., a);

  a = 12.01*g/mole;
  G4Element* elC  = new G4Element(name="Carbon",   symbol="C", z=6., a);

  a = 14.01*g/mole;
  G4Element* elN  = new G4Element(name="Nitrogen", symbol="N", z=7., a);

  a = 16.00*g/mole;
  G4Element* elO  = new G4Element(name="Oxygen",   symbol="O", z=8., a);

  a = 28.09*g/mole;
  G4Element* elSi = new G4Element(name="Silicon", symbol="Si", z=14., a);


//	------------- Materials -------------


// Air
// ---
  density = 1.29e-03*g/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);


// Aluminum 
// ---------
  a = 26.98*g/mole;
  density = 2.7*g/cm3;
  G4Material* Al ;
  Al = new G4Material(name="Aluminum", z=13., a, density); 


// Quartz
// -------
//  density = 2.200*g/cm3; // fused quartz 
  density = 2.64*g/cm3;  // crystalline quartz (c.f. PDG) 
  G4Material *Quartz = new G4Material(name="Quartz",density, nel=2);
  Quartz->AddElement(elSi, 1) ;
  Quartz->AddElement(elO , 2) ;


// PMMA C5H8O2 ( Acrylic )
// -------------
   density = 1.19*g/cm3;
   G4Material* Acrylic = new G4Material(name="Acrylic", density, nel=3);
   Acrylic->AddElement(elC, 5);
   Acrylic->AddElement(elH, 8);
   Acrylic->AddElement(elO, 2);


/////////////////////////////////////////////
// Construct Material Properties Tables
/////////////////////////////////////////////

  const G4int NUMENTRIES = 32;

  // Energy bins
  G4double X_RINDEX[NUMENTRIES] =
            { 2.034E-9*GeV, 2.068E-9*GeV, 2.103E-9*GeV, 2.139E-9*GeV,
              2.177E-9*GeV, 2.216E-9*GeV, 2.256E-9*GeV, 2.298E-9*GeV,
              2.341E-9*GeV, 2.386E-9*GeV, 2.433E-9*GeV, 2.481E-9*GeV,
              2.532E-9*GeV, 2.585E-9*GeV, 2.640E-9*GeV, 2.697E-9*GeV,
              2.757E-9*GeV, 2.820E-9*GeV, 2.885E-9*GeV, 2.954E-9*GeV,
              3.026E-9*GeV, 3.102E-9*GeV, 3.181E-9*GeV, 3.265E-9*GeV,
              3.353E-9*GeV, 3.446E-9*GeV, 3.545E-9*GeV, 3.649E-9*GeV,
              3.760E-9*GeV, 3.877E-9*GeV, 4.002E-9*GeV, 4.136E-9*GeV } ;


  // Air
  G4double RINDEX_AIR[NUMENTRIES] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 } ;

// Air refractive index at 20 oC and 1 atm (from PDG) 
  for(G4int j=0 ; j<NUMENTRIES ; j++){
    RINDEX_AIR[j] = RINDEX_AIR[j] + 2.73*std::pow(10.0,-4) ; 
    }

  G4MaterialPropertiesTable *MPT_Air = new G4MaterialPropertiesTable();
  MPT_Air->AddProperty("RINDEX", X_RINDEX, RINDEX_AIR, NUMENTRIES);
  Air->SetMaterialPropertiesTable(MPT_Air);

//////////////////////////////////////////////////////////////////////////////////////
//           Photomultiplier (PMT) window       
// The refractive index is for lime glass; 
// wavelength dependence is not included and value at 400nm is used.
//////////////////////////////////////////////////////////////////////////////////////

  // Refractive index 

  const G4int N_RINDEX_QUARTZ = 2 ;
  G4double X_RINDEX_QUARTZ[N_RINDEX_QUARTZ] = {0.0*eV, 10.0*eV};
  G4double RINDEX_QUARTZ[N_RINDEX_QUARTZ] = {1.54, 1.54};

  G4MaterialPropertiesTable *MPT_PMT = new G4MaterialPropertiesTable();
  MPT_PMT->AddProperty("RINDEX", X_RINDEX_QUARTZ, RINDEX_QUARTZ, N_RINDEX_QUARTZ);

  Quartz->SetMaterialPropertiesTable(MPT_PMT);


//////////////////////////////////////////////////////////////////
//               ACRYLIC Optical properties
//////////////////////////////////////////////////////////////////

// Refractive index 

  const G4int    N_RINDEX_ACRYLIC = 3 ;
  G4double X_RINDEX_ACRYLIC[N_RINDEX_ACRYLIC] = {320.0, 400.0, 500.0};  // Wavelength in nanometers
  G4double RINDEX_ACRYLIC[N_RINDEX_ACRYLIC] = {1.526, 1.507, 1.497};

    // Convert from nm to GeV

     for(G4int i=0;i<N_RINDEX_ACRYLIC; i++){
      X_RINDEX_ACRYLIC[i] = ((1239.84/X_RINDEX_ACRYLIC[i])*1E-9)*GeV;
      }

  G4MaterialPropertiesTable *MPT_Acrylic = new G4MaterialPropertiesTable();
  MPT_Acrylic->AddProperty("RINDEX", X_RINDEX_ACRYLIC, RINDEX_ACRYLIC, N_RINDEX_ACRYLIC);
  Acrylic->SetMaterialPropertiesTable(MPT_Acrylic);


//////////////////////////////////////////////////////////////////

  G4cout << *(G4Material::GetMaterialTable()) << G4endl ;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* UltraDetectorConstruction::ConstructMirror(G4VPhysicalVolume *World_phys){

  G4double Mirror_x = 40.0*cm;
  G4double Mirror_y = 40.0*cm;
  G4double Mirror_z = 1*cm;

  G4Box * boxMirror = new G4Box("Mirror",Mirror_x,Mirror_y,Mirror_z);

  // Get Air pointer from static funcion - (G4Material::GetMaterial)

G4String name;
G4Material *Al = G4Material::GetMaterial(name = "Aluminum");
G4LogicalVolume *logMirror ;
logMirror  = new G4LogicalVolume(boxMirror,Al,"Mirror",0,0,0);


G4ThreeVector SurfacePosition = G4ThreeVector(0*m,0*m,1.5*m) ;

// Rotate reflecting surface by 45. degrees around the OX axis.

G4RotationMatrix *Surfrot = new G4RotationMatrix(G4ThreeVector(1.0,0.0,0.0),-pi/4.);

G4VPhysicalVolume *physMirror ;
physMirror = new G4PVPlacement(Surfrot,SurfacePosition,"MirrorPV",logMirror,World_phys,false,0);

G4VisAttributes* SurfaceVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
SurfaceVisAtt->SetVisibility(true);
SurfaceVisAtt->SetForceWireframe(true);
logMirror->SetVisAttributes(SurfaceVisAtt);


//////////////////////////////////////////////////////////////////////////////////////////
//   Optical properties of the interface between the Air and Reflective Surface
//   For Mirror, reflectivity is set at 95% and specular reflection is assumed.


G4OpticalSurface *OpticalAirMirror = new G4OpticalSurface("AirMirrorSurface");
OpticalAirMirror->SetModel(unified);
OpticalAirMirror->SetType(dielectric_dielectric);
OpticalAirMirror->SetFinish(polishedfrontpainted);

const G4int NUM = 2;
G4double XX[NUM] = { 0.1E-9*GeV, 10.0E-9*GeV };
G4double ICEREFLECTIVITY[NUM]      = { 0.95, 0.95 };

G4MaterialPropertiesTable *AirMirrorMPT = new G4MaterialPropertiesTable();
AirMirrorMPT->AddProperty("REFLECTIVITY", XX, ICEREFLECTIVITY,NUM);
OpticalAirMirror->SetMaterialPropertiesTable(AirMirrorMPT);


G4LogicalBorderSurface *AirMirror ;
AirMirror = new G4LogicalBorderSurface("Air/Mirror Surface",World_phys,physMirror,OpticalAirMirror);

 return physMirror  ; 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* UltraDetectorConstruction::ConstructGround(G4VPhysicalVolume *World_phys){

  G4double Ground_x = 40.0*cm;
  G4double Ground_y = 40.0*cm;
  G4double Ground_z = 1*cm;

  G4Box * boxGround = new G4Box("Ground",Ground_x,Ground_y,Ground_z);

  // Get Air pointer from static funcion - (G4Material::GetMaterial)

G4String name;
G4Material *Al = G4Material::GetMaterial(name = "Aluminum");
G4LogicalVolume *logGround ;
logGround  = new G4LogicalVolume(boxGround,Al,"Ground",0,0,0);


G4ThreeVector SurfacePosition = G4ThreeVector(0*m,0*m,1.5*m) ;

// Rotate reflecting surface by 45. degrees around the OX axis.

G4RotationMatrix *Surfrot = new G4RotationMatrix(G4ThreeVector(1.0,0.0,0.0),-pi/4.);

G4VPhysicalVolume *physGround ;
physGround = new G4PVPlacement(Surfrot,SurfacePosition,"GroundPV",logGround,World_phys,false,0);

G4VisAttributes* SurfaceVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
SurfaceVisAtt->SetVisibility(true);
SurfaceVisAtt->SetForceWireframe(true);
logGround->SetVisAttributes(SurfaceVisAtt);


//////////////////////////////////////////////////////////////////////////////////////////
//   Optical properties of the interface between the Air and Reflective Surface
//   For Ground, reflectivity is set to 95% and diffusive reflection is assumed.


G4OpticalSurface *OpticalAirGround = new G4OpticalSurface("AirGroundSurface");
OpticalAirGround->SetModel(unified);
OpticalAirGround->SetType(dielectric_dielectric);
OpticalAirGround->SetFinish(groundfrontpainted);

 const G4int NUM = 2;
G4double XX[NUM] = { 0.1E-9*GeV, 10.0E-9*GeV };
G4double ICEREFLECTIVITY[NUM]      = { 0.95, 0.95 };

G4MaterialPropertiesTable *AirGroundMPT = new G4MaterialPropertiesTable();
AirGroundMPT->AddProperty("REFLECTIVITY", XX, ICEREFLECTIVITY,NUM);
OpticalAirGround->SetMaterialPropertiesTable(AirGroundMPT);


G4LogicalBorderSurface *AirGround ;
AirGround = new G4LogicalBorderSurface("Air/Ground Surface",World_phys,physGround,OpticalAirGround);

 return physGround  ; 

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* UltraDetectorConstruction::ConstructUVscope(G4VPhysicalVolume *World_phys){

//	------------- Volumes --------------

////////////////////////////////////////////////////////////////////////////////////////////////////////

G4cout << "#                                                    #" << G4endl ;  
G4cout << "#           Building the Telescope    ...            #" << G4endl ;  
G4cout << "#                                                    #" << G4endl ;  

/////////////////////////////////////////////////////////////
// UVscope housing is a cylinder made of 1 mm thick aluminum
/////////////////////////////////////////////////////////////

G4double UVscopeHeight    = 1030.0*mm ;
G4double UVscopeDiameter  = 518.0*mm ;
G4double UVscopeThickness = 1.0*mm   ;
G4double UVscopeBaffle    = 514.0*mm ; 

G4double UVscopeInnerRadius = UVscopeDiameter/2.0-UVscopeThickness ;
G4double UVscopeOuterRadius = UVscopeDiameter/2.0 ; 

G4ThreeVector UVscopePosition = G4ThreeVector(0.0*m,0.0*m,-1.0*m) ;
G4String name;
G4Material* Al = G4Material::GetMaterial(name = "Aluminum");


G4Tubs *solidUVscope = 
 new G4Tubs("UVscopeSolid",UVscopeInnerRadius,UVscopeOuterRadius,UVscopeHeight/2.0,0.0,twopi) ;
G4LogicalVolume *logicUVscope =
 new G4LogicalVolume(solidUVscope,Al,"UVscopeLV",0,0,0);
G4VPhysicalVolume *physicalUVscope =
 new G4PVPlacement(0,UVscopePosition,"UVSCopePV",logicUVscope,World_phys,false,0);


//////////////////////////////////////
// Back cover of the UVscope cylinder
//////////////////////////////////////

G4Tubs *solidUVscopeBack = 
 new G4Tubs("UVscopeBackSolid",0.0,UVscopeOuterRadius,UVscopeThickness/2.0,0.0,twopi) ;

G4LogicalVolume *logicUVscopeBack = 
 new G4LogicalVolume(solidUVscopeBack,Al,"UVscopeBackLV",0,0,0);

G4ThreeVector UVscopeBackPosition ;
UVscopeBackPosition =  UVscopePosition+G4ThreeVector(0.0*mm,0.0*mm,-(UVscopeHeight/2.0+UVscopeThickness/2.0)) ;
G4VPhysicalVolume *physicalUVscopeBack = 
 new G4PVPlacement(0,UVscopeBackPosition,"UVscopeBack",logicUVscopeBack,World_phys,false,0);



////////////////////////////////////////////////////////////////////////////////////////////////////////

  G4cout << "#                                                    #" << G4endl ;  
  G4cout << "#           Building the Fresnel lens ...            #" << G4endl ;  
  G4cout << "#                                                    #" << G4endl ;  

G4double      LensDiameter        = 457*mm ; // Size of the optical active area of the lens.
G4int      LensNumOfGrooves    = 13 ;
//G4int      LensNumOfGrooves    = 129 ;
//G4int      LensNumOfGrooves    = 1287 ;

G4double      LensBorderThickness = 2.8*mm ;     // Thickness of the border area. 
G4double      LensFocalLength     = 441.973*mm ; // This parameter depends on the lens geometry, etc !!
G4Material   *LensMaterial        = G4Material::GetMaterial(name = "Acrylic") ;
G4ThreeVector LensPosition        = UVscopePosition+G4ThreeVector(0.0*mm,0.0*mm,UVscopeHeight/2.0-UVscopeBaffle) ;


UltraFresnelLens *FresnelLens = new UltraFresnelLens(LensDiameter,LensNumOfGrooves,LensMaterial,World_phys,LensPosition) ;


///////////////////////////////////
// Lens supporting ring (aluminum)
///////////////////////////////////

G4Tubs *solidLensFrame = new G4Tubs("LensFrame",LensDiameter/2.0,UVscopeInnerRadius,LensBorderThickness/2.0,0.0,twopi) ;
G4LogicalVolume *logicLensFrame = new G4LogicalVolume(solidLensFrame,Al,"LensFrameLV",0,0,0);

G4ThreeVector LensFramePosition ;
LensFramePosition = LensPosition+G4ThreeVector(0.0*mm,0.0*mm,-((FresnelLens->GetThickness())/2.0+solidLensFrame->GetDz())) ;

G4VPhysicalVolume *physicalLensFrame =
  new G4PVPlacement(0,LensFramePosition,"LensFramePV",logicLensFrame,World_phys,false,0);

////////////////////////////////////////////////////////////////////////////////////////////////////////


  G4cout << "#                                                    #" << G4endl ;  
  G4cout << "#         Building the photomultiplier ...           #" << G4endl ;  
  G4cout << "#                                                    #" << G4endl ;  


// Photomultiplier window is a spherical section made of quartz

G4double PMT_thick   =   1.0*mm ; // Thickness of PMT window
G4double PMT_curv    =  65.5*mm ; // Radius of curvature of PMT window
G4double StartTheta  = (180.0-31.2)*pi/180. ;
G4double EndTheta    = 31.2*pi/180. ;

G4Sphere *solidPMT ;
solidPMT = new G4Sphere("PMT_solid",PMT_curv-PMT_thick,PMT_curv,0.0,twopi,StartTheta,EndTheta);

G4Material* Quartz = G4Material::GetMaterial(name = "Quartz");
G4LogicalVolume * logicalPMT ;
logicalPMT = new G4LogicalVolume(solidPMT,Quartz,"PMT_log",0,0,0);


// Place PMT is at Lens Focus

G4ThreeVector PMTpos = LensPosition + G4ThreeVector(0.0*cm,0.0*cm,-(LensFocalLength+PMT_curv)) ;

// Rotate PMT window through the axis OX by an angle = 180. degrees

G4RotationMatrix *PMTrot = new G4RotationMatrix(G4ThreeVector(1.0,0.0,0.0),pi);
G4VPhysicalVolume *physPMT ;
physPMT  = new G4PVPlacement(PMTrot,PMTpos,"PMT1",logicalPMT,World_phys,false,0);

  if(!PMTSD)
    {
      PMTSD = new UltraPMTSD("PMTSD");
      SDmanager->AddNewDetector( PMTSD );
    }

  if (logicalPMT){logicalPMT->SetSensitiveDetector(PMTSD);}

G4VisAttributes* PMTVisAtt   = new G4VisAttributes(true,G4Colour(0.0,0.0,1.0)) ;   
logicalPMT->SetVisAttributes(PMTVisAtt);

//////////////////////////////////////////////////////////////////////////////////////////
//   Optical properties of the interface between the Air and the walls of the 
//   UVscope cylinder (5% reflectivity)


  G4cout << "#    Defining interface's optical properties  ...    #" << G4endl ;  
  G4cout << "#                                                    #" << G4endl ;  


G4OpticalSurface *OpticalAirPaint = new G4OpticalSurface("AirPaintSurface");
OpticalAirPaint->SetModel(unified);
OpticalAirPaint->SetType(dielectric_dielectric);
OpticalAirPaint->SetFinish(groundfrontpainted);

const G4int NUM = 2;
G4double XX[NUM] = { 2.030E-9*GeV, 4.144E-9*GeV };
G4double BLACKPAINTREFLECTIVITY[NUM]      = { 0.05, 0.05 };
//G4double WHITEPAINTREFLECTIVITY[NUM]      = { 0.99, 0.99 };

G4MaterialPropertiesTable *AirPaintMPT = new G4MaterialPropertiesTable();
AirPaintMPT->AddProperty("REFLECTIVITY", XX, BLACKPAINTREFLECTIVITY,NUM);
OpticalAirPaint->SetMaterialPropertiesTable(AirPaintMPT);

//OpticalAirPaint->DumpInfo();

G4LogicalBorderSurface *AirCylinder ;
AirCylinder = new G4LogicalBorderSurface("Air/UVscope Cylinder Surface",World_phys,physicalUVscope,OpticalAirPaint);

G4LogicalBorderSurface *AirLensFrame ;
AirLensFrame = new G4LogicalBorderSurface("Air/LensFrame Surface",World_phys,physicalLensFrame,OpticalAirPaint);

G4LogicalBorderSurface *AirBackCover ;
AirBackCover = new G4LogicalBorderSurface("Air/UVscope Back Cover Surface",World_phys,physicalUVscopeBack,OpticalAirPaint);


/////////////////////////////////////////////////////////////////////////////////////


   G4VisAttributes* LensVisAtt  = new G4VisAttributes(G4Colour(1.0,0.0,0.0)) ;   // Red
   LensVisAtt ->SetVisibility(true);


   if (FresnelLens){
   FresnelLens->GetPhysicalVolume()->GetLogicalVolume()->SetVisAttributes(LensVisAtt);
   }

   G4VisAttributes* UVscopeVisAtt  = new G4VisAttributes(G4Colour(0.5,0.5,0.5)) ;   // Gray
   UVscopeVisAtt ->SetVisibility(true);

   physicalUVscope     ->GetLogicalVolume()->SetVisAttributes(UVscopeVisAtt);
   physicalUVscopeBack ->GetLogicalVolume()->SetVisAttributes(UVscopeVisAtt);
   physicalLensFrame   ->GetLogicalVolume()->SetVisAttributes(UVscopeVisAtt);

/////////////////////////////////////////////////////////////////////////////////////

  G4cout << "#                                                    #" << G4endl ;  
  G4cout << "#               UVscope is built ! ...               #" << G4endl ;  
  G4cout << "#                                                    #" << G4endl ;  

  return physicalUVscope;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



