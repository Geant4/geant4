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
#include "UltraDetectorMessenger.hh"
#include "UltraPMTSD.hh"
#include "UltraFresnelLens.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4GeometryManager.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
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
#include "G4Log.hh"
#include "G4SDManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UltraDetectorConstruction::UltraDetectorConstruction() :
  fReflectorOpticalSurface(0),
  logicalPMT(0),
  fReflectorLog(0),
  fIsReflectorConstructed(false)
{
  // Define wavelength limits for materials definition
  lambda_min = 200*nm ;
  lambda_max = 700*nm ;

  fDetectorMessenger = new UltraDetectorMessenger(this);

  ConstructTableMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UltraDetectorConstruction::~UltraDetectorConstruction()
{
  delete fDetectorMessenger;

  delete fReflectorOpticalSurface;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* UltraDetectorConstruction::Construct()
{
//	The experimental Hall
//	---------------------

  auto World_x = 1.*m;
  auto World_y = 1.*m;
  auto World_z = 2*m;

  auto World_box = new G4Box("World",World_x,World_y,World_z);

  // Get Air pointer from static funcion - (G4Material::GetMaterial)
  auto Air = G4Material::GetMaterial("Air");
  auto World_log = new G4LogicalVolume(World_box,Air,"World",0,0,0);

  fWorld_phys   = new G4PVPlacement(0,G4ThreeVector(),"World",World_log,0,false,0);

   auto UniverseVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
   UniverseVisAtt->SetVisibility(true);
   UniverseVisAtt->SetForceWireframe(true);
   World_log->SetVisAttributes(UniverseVisAtt);
   World_log->SetVisAttributes (G4VisAttributes::GetInvisible());



  G4cout << "\n \n \n \n \n \n \n \n \n \n \n \n \n " << G4endl ;

  G4cout << "######################################################" << G4endl ;
  G4cout << "#                                                    #" << G4endl ;
  G4cout << "#                                                    #" << G4endl ;
  G4cout << "#          UltraDetectorConstruction:                #" << G4endl ;
  G4cout << "#                                                    #" << G4endl ;
  G4cout << "#                                                    #" << G4endl ;

  ConstructUVscope();


  G4cout << "#                                                    #" << G4endl ;
  G4cout << "#                                                    #" << G4endl ;
  G4cout << "######################################################" << G4endl ;

  fIsReflectorConstructed = false;

  return fWorld_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UltraDetectorConstruction::ConstructSDandField()
{
  auto PMTSD = new UltraPMTSD("PMTSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(PMTSD);
  SetSensitiveDetector(logicalPMT,PMTSD);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UltraDetectorConstruction::ConstructTableMaterials()
{
  G4double a, z, density;

//	------------- Elements -------------
  a = 1.01*g/mole;
  auto elH  = new G4Element("Hydrogen", "H", z=1., a);

  a = 12.01*g/mole;
  auto elC  = new G4Element("Carbon", "C", z=6., a);

  a = 14.01*g/mole;
  auto elN  = new G4Element("Nitrogen", "N", z=7., a);

  a = 16.00*g/mole;
  auto elO  = new G4Element("Oxygen",  "O", z=8., a);

  a = 28.09*g/mole;
  auto elSi = new G4Element("Silicon", "Si", z=14., a);


//	------------- Materials -------------


// Air
// ---
  density = 1.29e-03*g/cm3;
  auto Air = new G4Material("Air", density, 2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);


// Aluminum
// ---------
  a = 26.98*g/mole;
  density = 2.7*g/cm3;
  new G4Material("Aluminum", z=13., a, density);


// Quartz
// -------
//  density = 2.200*g/cm3; // fused quartz
  density = 2.64*g/cm3;  // crystalline quartz (c.f. PDG)
  auto Quartz = new G4Material("Quartz",density, 2);
  Quartz->AddElement(elSi, 1) ;
  Quartz->AddElement(elO , 2) ;


// PMMA C5H8O2 ( Acrylic )
// -------------
   density = 1.19*g/cm3;
   auto Acrylic = new G4Material("Acrylic", density, 3);
   Acrylic->AddElement(elC, 5);
   Acrylic->AddElement(elH, 8);
   Acrylic->AddElement(elO, 2);


/////////////////////////////////////////////
// Construct Material Properties Tables
/////////////////////////////////////////////

  // Energy bins
  std::vector<G4double> X_RINDEX = {h_Planck*c_light/lambda_max, h_Planck*c_light/lambda_min} ;


  // Air
  std::vector<G4double> RINDEX_AIR{1.00, 1.00} ;
// Air refractive index at 20 oC and 1 atm (from PDG)
  for(auto&& i : RINDEX_AIR){
    i = i  + 2.73*std::pow(10.0,-4) ;
  }

  auto MPT_Air = new G4MaterialPropertiesTable();
  MPT_Air->AddProperty("RINDEX", X_RINDEX, RINDEX_AIR);
  Air->SetMaterialPropertiesTable(MPT_Air);

//////////////////////////////////////////////////////////////////////////////////////
//           Photomultiplier (PMT) window
// The refractive index is for lime glass;
// wavelength dependence is not included and value at 400nm is used.
//////////////////////////////////////////////////////////////////////////////////////

  // Refractive index

  std::vector<G4double> X_RINDEX_QUARTZ{h_Planck*c_light/lambda_max, h_Planck*c_light/lambda_min} ;
  std::vector<G4double> RINDEX_QUARTZ{1.54, 1.54};

  auto MPT_PMT = new G4MaterialPropertiesTable();
  MPT_PMT->AddProperty("RINDEX", X_RINDEX_QUARTZ, RINDEX_QUARTZ);

  Quartz->SetMaterialPropertiesTable(MPT_PMT);


//////////////////////////////////////////////////////////////////
//               ACRYLIC Optical properties
//////////////////////////////////////////////////////////////////

// Refractive index

  const auto NENTRIES = 11 ;

  std::vector<G4double> RINDEX_ACRYLIC;
  std::vector<G4double> ENERGY_ACRYLIC;

// Parameterization for refractive index of High Grade PMMA

  G4double bParam[4] = {1760.7010,-1.3687,2.4388e-3,-1.5178e-6} ;
  auto lambda = 0.;

  for(auto i=0;i<NENTRIES; ++i){

    // want energies in increasing order
    lambda = lambda_min + (NENTRIES-1-i) * (lambda_max-lambda_min)/float(NENTRIES-1);
    RINDEX_ACRYLIC.push_back(0.0);

    for (auto jj=0 ; jj<4 ; jj++)
    {
      RINDEX_ACRYLIC[i] +=  (bParam[jj]/1000.0)*std::pow(lambda/nm,jj) ;
    }

    ENERGY_ACRYLIC.push_back(h_Planck*c_light/lambda);  // Convert from wavelength to energy ;
//  G4cout << ENERGY_ACRYLIC[i]/eV << " " << lambda/nm << " " << RINDEX_ACRYLIC[i] << G4endl ;

  }

  auto MPT_Acrylic = new G4MaterialPropertiesTable();
  MPT_Acrylic->AddProperty("RINDEX", ENERGY_ACRYLIC, RINDEX_ACRYLIC);

// Absorption
  std::vector<G4double> LAMBDAABS
  {
    100.0,
    246.528671, 260.605103, 263.853516, 266.019104, 268.726105,
    271.433136, 273.598724, 276.305725, 279.554138, 300.127380,
    320.159241, 340.191101, 360.764343, 381.337585, 399.745239,
    421.401276, 440.891724, 460.382172, 480.414001, 500.987274,
    520.477722, 540.509583, 559.458618,
    700.0
  } ;

  std::vector<G4double> ABS   // Transmission (in %) of  3mm thick PMMA
  {
    0.0000000,
    0.0000000,  5.295952,  9.657321, 19.937695, 29.283491,
    39.252335, 48.598133, 58.255451, 65.109039, 79.439247,
    85.669785, 89.719627, 91.277260, 91.588783, 91.900307,
    91.588783, 91.277260, 91.277260, 91.588783, 91.588783,
    91.900307, 91.900307, 91.588783,
    91.5
  } ;


  MPT_Acrylic->AddProperty("ABSLENGTH", new G4MaterialPropertyVector()) ;
  for(size_t i=0;i<ABS.size(); ++i){
    auto energy    = h_Planck*c_light/(LAMBDAABS[i]*nm) ;
    auto abslength = 0.0;

    if (ABS[i] <= 0.0) {
      abslength = 1.0/kInfinity ;
    }
    else {
      abslength = -3.0*mm/(G4Log(ABS[i]/100.0)) ;
    }

    MPT_Acrylic->AddEntry("ABSLENGTH", energy, abslength);

  }

  Acrylic->SetMaterialPropertiesTable(MPT_Acrylic);


//////////////////////////////////////////////////////////////////

  G4cout << *(G4Material::GetMaterialTable()) << G4endl ;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UltraDetectorConstruction::ConstructReflector()
{
  const auto x = 40.0*cm;
  const auto y = 40.0*cm;
  const auto z = 1*cm;

  auto box = new G4Box("Mirror",x,y,z);

  // Get Air pointer from static funcion - (G4Material::GetMaterial)
  auto Al = G4Material::GetMaterial("Aluminum");

  fReflectorLog = new G4LogicalVolume(box,Al,"Reflector",0,0,0);

  auto SurfacePosition = G4ThreeVector(0*m,0*m,1.5*m) ;

  // Rotate reflecting surface by 45. degrees around the OX axis.

  auto Surfrot = new G4RotationMatrix(G4ThreeVector(1.0,0.0,0.0),-pi/4.);

  new G4PVPlacement(Surfrot,SurfacePosition,"MirrorPV",fReflectorLog,fWorld_phys,false,0);

  auto SurfaceVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  SurfaceVisAtt->SetVisibility(true);
  SurfaceVisAtt->SetForceWireframe(true);
  fReflectorLog->SetVisAttributes(SurfaceVisAtt);

  fReflectorOpticalSurface = new G4OpticalSurface("ReflectorOpticalSurface");
  fReflectorOpticalSurface->SetModel(unified);
  fReflectorOpticalSurface->SetType(dielectric_dielectric);

  std::vector<G4double> XX{h_Planck*c_light/lambda_max, h_Planck*c_light/lambda_min} ;
  std::vector<G4double> ICEREFLECTIVITY{ 0.95, 0.95 };

  auto AirMirrorMPT = new G4MaterialPropertiesTable();
  AirMirrorMPT->AddProperty("REFLECTIVITY", XX, ICEREFLECTIVITY);
  fReflectorOpticalSurface->SetMaterialPropertiesTable(AirMirrorMPT);

  new G4LogicalSkinSurface("ReflectorSurface",fReflectorLog,fReflectorOpticalSurface);

#ifdef G4MULTITHREADED
  auto runManager = G4MTRunManager::GetMasterRunManager();
  //runManager->SetNumberOfThreads(2);
#else
  auto runManager = G4RunManager::GetRunManager();
#endif

  runManager->GeometryHasBeenModified();

  fIsReflectorConstructed = true;
}


void UltraDetectorConstruction::SetReflectorOpticalProperties()
{
  if (fReflectionType == "ground") {
    G4cout << "Using ground reflecting surface " << G4endl ;
    if (fReflectorOpticalSurface)
      fReflectorOpticalSurface->SetFinish(groundfrontpainted);
  }
  else {
    G4cout << "Using mirror reflecting surface " << G4endl ;
    if (fReflectorOpticalSurface)
      fReflectorOpticalSurface->SetFinish(polishedfrontpainted);
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UltraDetectorConstruction::ConstructUVscope()
{

  //	------------- Volumes --------------

  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  G4cout << "#                                                    #" << G4endl ;
  G4cout << "#           Building the Telescope    ...            #" << G4endl ;
  G4cout << "#                                                    #" << G4endl ;

  /////////////////////////////////////////////////////////////
  // UVscope housing is a cylinder made of 1 mm thick aluminum
  /////////////////////////////////////////////////////////////

  auto UVscopeHeight    = 1030.0*mm ;
  auto UVscopeDiameter  = 518.0*mm ;
  auto UVscopeThickness = 1.0*mm   ;
  auto UVscopeBaffle    = 514.0*mm ;

  auto UVscopeInnerRadius = UVscopeDiameter/2.0-UVscopeThickness ;
  auto UVscopeOuterRadius = UVscopeDiameter/2.0 ;

  auto UVscopePosition = G4ThreeVector(0.0*m,0.0*m,-1.0*m) ;
  auto Al = G4Material::GetMaterial("Aluminum");

  auto solidUVscope =
    new G4Tubs("UVscopeSolid",UVscopeInnerRadius,UVscopeOuterRadius,UVscopeHeight/2.0,0.0,twopi) ;
  auto logicUVscope =
    new G4LogicalVolume(solidUVscope,Al,"UVscopeLV",0,0,0);
  auto physicalUVscope =
    new G4PVPlacement(0,UVscopePosition,"UVSCopePV",logicUVscope,fWorld_phys,false,0);


  //////////////////////////////////////
  // Back cover of the UVscope cylinder
  //////////////////////////////////////

  auto solidUVscopeBack =
    new G4Tubs("UVscopeBackSolid",0.0,UVscopeOuterRadius,UVscopeThickness/2.0,0.0,twopi) ;

  auto logicUVscopeBack =
    new G4LogicalVolume(solidUVscopeBack,Al,"UVscopeBackLV",0,0,0);

  auto UVscopeBackPosition  =
    UVscopePosition+G4ThreeVector(0.0*mm,0.0*mm,-(UVscopeHeight/2.0+UVscopeThickness/2.0)) ;
  auto physicalUVscopeBack =
    new G4PVPlacement(0,UVscopeBackPosition,"UVscopeBack",logicUVscopeBack,fWorld_phys,false,0);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  G4cout << "#                                                    #" << G4endl ;
  G4cout << "#           Building the Fresnel lens ...            #" << G4endl ;
  G4cout << "#                                                    #" << G4endl ;

  auto      LensDiameter        = 457*mm ; // Size of the optical active area of the lens.
  auto      LensNumOfGrooves    = 13 ;
  //auto      LensNumOfGrooves    = 129 ;
  //auto      LensNumOfGrooves    = 1287 ;

  auto      LensBorderThickness = 2.8*mm ;     // Thickness of the border area.
  auto      LensFocalLength     = 441.973*mm ; // This parameter depends on the lens geometry, etc !!
  auto      LensMaterial        = G4Material::GetMaterial("Acrylic") ;
  auto      LensPosition        = UVscopePosition+G4ThreeVector(0.0*mm,0.0*mm,UVscopeHeight/2.0-UVscopeBaffle) ;


  FresnelLens = new UltraFresnelLens(LensDiameter,LensNumOfGrooves,LensMaterial,fWorld_phys,LensPosition) ;


  ///////////////////////////////////
  // Lens supporting ring (aluminum)
  ///////////////////////////////////

  auto solidLensFrame = new G4Tubs("LensFrame",LensDiameter/2.0,UVscopeInnerRadius,LensBorderThickness/2.0,0.0,twopi) ;
  auto logicLensFrame = new G4LogicalVolume(solidLensFrame,Al,"LensFrameLV",0,0,0);

  auto LensFramePosition = LensPosition+G4ThreeVector(0.0*mm,0.0*mm,-((FresnelLens->GetThickness())/2.0+solidLensFrame->GetZHalfLength())) ;

  auto physicalLensFrame =
    new G4PVPlacement(0,LensFramePosition,"LensFramePV",logicLensFrame,fWorld_phys,false,0);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////


  G4cout << "#                                                    #" << G4endl ;
  G4cout << "#         Building the photomultiplier ...           #" << G4endl ;
  G4cout << "#                                                    #" << G4endl ;


  // Photomultiplier window is a spherical section made of quartz

  auto PMT_thick   =   1.0*mm ; // Thickness of PMT window
  auto PMT_curv    =  65.5*mm ; // Radius of curvature of PMT window
  auto StartTheta  = (180.0-31.2)*pi/180. ;
  auto EndTheta    = 31.2*pi/180. ;

  auto solidPMT =
    new G4Sphere("PMT_solid",PMT_curv-PMT_thick,PMT_curv,0.0,twopi,StartTheta,EndTheta);

  auto Quartz = G4Material::GetMaterial("Quartz");
  logicalPMT = new G4LogicalVolume(solidPMT,Quartz,"PMT_log",0,0,0);


  // Place PMT is at Lens Focus

  auto PMTpos = LensPosition + G4ThreeVector(0.0*cm,0.0*cm,-(LensFocalLength+PMT_curv)) ;

  // Rotate PMT window through the axis OX by an angle = 180. degrees

  auto PMTrot = new G4RotationMatrix(G4ThreeVector(1.0,0.0,0.0),pi);
  new G4PVPlacement(PMTrot,PMTpos,"PMT1",logicalPMT,fWorld_phys,false,0);


  auto PMTVisAtt   = new G4VisAttributes(true,G4Colour(0.0,0.0,1.0)) ;
  logicalPMT->SetVisAttributes(PMTVisAtt);

  //////////////////////////////////////////////////////////////////////////////////////////
  //   Optical properties of the interface between the Air and the walls of the
  //   UVscope cylinder (5% reflectivity)


  G4cout << "#    Defining interface's optical properties  ...    #" << G4endl ;
  G4cout << "#                                                    #" << G4endl ;


  auto OpticalAirPaint = new G4OpticalSurface("AirPaintSurface");
  OpticalAirPaint->SetModel(unified);
  OpticalAirPaint->SetType(dielectric_dielectric);
  OpticalAirPaint->SetFinish(groundfrontpainted);

  std::vector<G4double> XX = {h_Planck*c_light/lambda_max, h_Planck*c_light/lambda_min} ;
  std::vector<G4double> BLACKPAINTREFLECTIVITY      = { 0.05, 0.05 };
  //std::vector<G4double> WHITEPAINTREFLECTIVITY      = { 0.99, 0.99 };

  auto AirPaintMPT = new G4MaterialPropertiesTable();
  AirPaintMPT->AddProperty("REFLECTIVITY", XX, BLACKPAINTREFLECTIVITY);
  OpticalAirPaint->SetMaterialPropertiesTable(AirPaintMPT);

  //OpticalAirPaint->DumpInfo();

  new G4LogicalBorderSurface("Air/UVscope Cylinder Surface",fWorld_phys,physicalUVscope,OpticalAirPaint);

  new G4LogicalBorderSurface("Air/LensFrame Surface",fWorld_phys,physicalLensFrame,OpticalAirPaint);

  new G4LogicalBorderSurface("Air/UVscope Back Cover Surface",fWorld_phys,physicalUVscopeBack,OpticalAirPaint);


  /////////////////////////////////////////////////////////////////////////////////////

  auto LensVisAtt  = new G4VisAttributes(G4Colour(1.0,0.0,0.0)) ;   // Red
  LensVisAtt ->SetVisibility(true);


  if (FresnelLens){
    FresnelLens->GetPhysicalVolume()->GetLogicalVolume()->SetVisAttributes(LensVisAtt);
  }

  auto UVscopeVisAtt  = new G4VisAttributes(G4Colour(0.5,0.5,0.5)) ;   // Gray
  UVscopeVisAtt ->SetVisibility(true);

  physicalUVscope     ->GetLogicalVolume()->SetVisAttributes(UVscopeVisAtt);
  physicalUVscopeBack ->GetLogicalVolume()->SetVisAttributes(UVscopeVisAtt);
  physicalLensFrame   ->GetLogicalVolume()->SetVisAttributes(UVscopeVisAtt);

  /////////////////////////////////////////////////////////////////////////////////////

  G4cout << "#                                                    #" << G4endl ;
  G4cout << "#               UVscope is built ! ...               #" << G4endl ;
  G4cout << "#                                                    #" << G4endl ;

}


void UltraDetectorConstruction::SetReflectionType(G4String rtype)
{
#ifdef G4MULTITHREADED
  auto runManager = G4MTRunManager::GetMasterRunManager();
  //runManager->SetNumberOfThreads(2);
#else
  auto runManager = G4RunManager::GetRunManager();
#endif

  fReflectionType = rtype;

  if (fReflectionType == "none") {
    if (fIsReflectorConstructed) {
      // Cleanup old geometry to delete reflecting surface
      runManager->ReinitializeGeometry(true);
    }
  }
  else {
    if (!fIsReflectorConstructed) {
      ConstructReflector();
    }
    SetReflectorOpticalProperties();
  }
}
