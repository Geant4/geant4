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
// Author: F. Poignant, floriane.poignant@gmail.com
//
//
//
//
//    ******************************************
//    *                                        *
//    *    STCyclotronDetectorConstruction.cc  *
//    *                                        *
//    ******************************************
//
//
#include "STCyclotronDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4NistManager.hh" 

#include "G4Element.hh" 
#include "G4Material.hh"

#include "G4Box.hh"
#include "G4Tubs.hh" 
#include "G4Cons.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Region.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "STCyclotronRun.hh"
#include "STCyclotronSensitiveTarget.hh"
#include "STCyclotronSensitiveFoil.hh"
#include "STCyclotronDetectorMessenger.hh"
#include "STCyclotronRunAction.hh"


STCyclotronDetectorConstruction::STCyclotronDetectorConstruction()
 :fTarget_diameter(0),fIsotopeName(0),fIsotopeZ(0),fIsotopeN(0),fIsotopeA(0), 
  fElementName(0),fElementSymbole(0),fElementNComponents(0),fElementAbundance(0),
  fNaturalElementName(0),fNaturalMaterialFractionMass(0), fDensity_target(0), 
  fTarget_NComponents(0), fMaterialFractionMass(0),fIsotopeNameFoil(0),
  fIsotopeZFoil(0),fIsotopeNFoil(0),fIsotopeAFoil(0), fElementNameFoil(0),
  fElementSymboleFoil(0),fElementNComponentsFoil(0),fElementAbundanceFoil(0),
  fNaturalElementNameFoil(0),fNaturalMaterialFractionMassFoil(0), 
  fDensity_foil(0), fFoil_NComponents(0), fMaterialFractionMassFoil(0), 
  fTarget_thickness(0), fFoil_thickness(0), fTarget_Material(0), fFoil_Material(0),
  fZ_foil_position(0), fSolidFoil(nullptr),fLogicFoil(nullptr), fPhysFoil(nullptr),
  fLogicWorld(nullptr),
  fLayer_z_position_PART3(0), fPhysLayer_PART3(nullptr), fPhysTube_PART3(nullptr),
  fTube_outerRadius_PART4(0), fTube_length_PART4(0),fLayer_z_position_PART4(0), 
  fPhysTube_PART4(nullptr), fPhysLayer_PART4(nullptr), 
  fLayer1_z_position_PART4(0),fPhysLayer1_PART4(nullptr),
  fLogicTarget(nullptr), fTarget_z_position(0),fSolidTarget(nullptr), 
  fPhysTarget(nullptr),
  fLayer1_z_position_PART5(0), fPhysLayer1_PART5(nullptr), 
  fLayer2_z_position_PART5(0), fPhysLayer2_PART5(nullptr), 
  fLayer3_z_position_PART5(0),fPhysLayer3_PART5(nullptr),  
  fRegionTarget(nullptr), fRegionFoil(nullptr), fTargetVolume(0), fFoilVolume(0)
{ 
  fDetectorMessenger = new STCyclotronDetectorMessenger(this);
}

STCyclotronDetectorConstruction::~STCyclotronDetectorConstruction()
{
  delete fDetectorMessenger;
}

G4VPhysicalVolume* STCyclotronDetectorConstruction::Construct()
{  
  //Initialization of messenger parameters
  fTarget_diameter = 7.*mm;
  fDensity_target = 8.9*g/cm3;
  fTarget_thickness = 0.35*mm;
  fFoil_thickness = 0.000001*mm;
  fDensity_foil = 2.7*g/cm3;
   

  //Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;


  //Create the world
  G4double world_hx = 1.*m;
  G4double world_hy = 1.*m;
  G4double world_hz = 1.*m;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  G4Box* solidWorld 
    = new G4Box("World",
		world_hx, 
		world_hy, 
		world_hz);
  
  fLogicWorld
    = new G4LogicalVolume(solidWorld,
			  world_mat,
			  "World");

  G4VPhysicalVolume* physWorld 
    = new G4PVPlacement(0,                     //no rotation
			G4ThreeVector(),       //at (0,0,0)
			fLogicWorld,            //its logical volume
			"World",               //its name
			0,                     //its mother  volume
			false,                 //no boolean operation
			0);                    //copy number



  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////////Create the detector////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  

  //Overall parameters
  G4double startAngle     = 0.*deg;
  G4double spanningAngle  = 360.*deg;
  
  ////////////////////////////////////
  /////////Define materials///////////
  ////////////////////////////////////


  //ALUMINIUM//
 
  G4Material* al = nist->FindOrBuildMaterial("G4_Al");



  //Create vacuum around the beam//
  G4double vacuum_atomic_number, vacuum_mass_of_mole, vacuum_density,vacuum_pressure,vacuum_temperature;
  vacuum_atomic_number = 1.;
  vacuum_mass_of_mole = 1.008*g/mole;
  vacuum_density      = 1.e-30*g/cm3;
  vacuum_pressure     = 1.e-8*bar;
  vacuum_temperature  = 293.*kelvin;     //from PhysicalConstants.h
  G4Material* vacuum_beam = new G4Material("vacuumBeam", vacuum_atomic_number,
					   vacuum_mass_of_mole, vacuum_density,
					   kStateGas,vacuum_temperature,
					   vacuum_pressure);



  //HELIUM
  G4double helium_Z, helium_A, helium_density, helium_pressure, helium_temperature;
  helium_Z = 2.;
  helium_A = 4*g/mole;
  helium_density = 0.1785e-3*g/cm3; // with T=0°, 1 atm. To modify !
  helium_pressure = 2.*bar;
  helium_temperature = 293.*kelvin; //15 to 20°
  G4Material* helium = new G4Material("helium", helium_Z, helium_A, helium_density, kStateGas, helium_temperature, helium_pressure);


  //PLATINIUM
  G4Material* pt =  nist->FindOrBuildMaterial("G4_Pt");


  //TARGET MATERIAL INITIALIZATION
  /*G4String name;
  G4String symbole;
  G4int ncomponents;
  G4int n;
  G4double z_isotope;
  G4double abundance;
  G4double fractionmass;
  G4double a;*/
 

  /*  
  //Pure Ni64
  
  G4Isotope* Ni64 = new G4Isotope(name="Zi64", z_isotope=28., n=64, a=64.*g/mole);
  
  G4Element* pureNi64 = new G4Element(name="pureNi64",symbole="64Ni", ncomponents=1);
  pureNi64->AddIsotope(Ni64, abundance = 100.*perCent);


  fTarget_Material = new G4Material("FTarget_Material",fDensity_target,ncomponents=1);
  fTarget_Material->AddElement(pureNi64, fractionmass=1.);
 */
  /*  
  //Ni64 94% enriched - Sz

  G4Isotope* Ni64 = new G4Isotope(name="Zi64", z_isotope=28., n=64, a=64.*g/mole);
  G4Isotope* Ni58 = new G4Isotope(name="Zi58", z_isotope=28., n=58, a=58.*g/mole);
  G4Isotope* Ni60 = new G4Isotope(name="Zi60", z_isotope=28., n=60, a=60.*g/mole);
  G4Isotope* Ni61 = new G4Isotope(name="Zi61", z_isotope=28., n=61, a=61.*g/mole);
  G4Isotope* Ni62 = new G4Isotope(name="Zi62", z_isotope=28., n=62, a=62.*g/mole);


  G4Element* Ni64enriched95 = new G4Element(name="Ni64enriched95",symbole="64Ni_95", ncomponents=5);
  Ni64enriched95->AddIsotope(Ni64, abundance = 95.*perCent);
  Ni64enriched95->AddIsotope(Ni58, abundance = 2.6*perCent);
  Ni64enriched95->AddIsotope(Ni60, abundance = 1.72*perCent);
  Ni64enriched95->AddIsotope(Ni61, abundance = 0.15*perCent);
  Ni64enriched95->AddIsotope(Ni62, abundance = 0.53*perCent);


  fTarget_Material = new G4Material("FTarget_Material",fDensity_target,ncomponents=1);
  fTarget_Material->AddElement(Ni64enriched95, fractionmass=1.);

  //fTarget_Material =  nist->FindOrBuildMaterial("G4_Y");
  */
  
  /*

//Ni64 95% enriched - Obata

  G4Isotope* Ni64 = new G4Isotope(name="Zi64", z_isotope=28., n=64, a=64.*g/mole);
  G4Isotope* Ni58 = new G4Isotope(name="Zi58", z_isotope=28., n=58, a=58.*g/mole);
  G4Isotope* Ni60 = new G4Isotope(name="Zi60", z_isotope=28., n=60, a=60.*g/mole);
  G4Isotope* Ni61 = new G4Isotope(name="Zi61", z_isotope=28., n=61, a=61.*g/mole);
  G4Isotope* Ni62 = new G4Isotope(name="Zi62", z_isotope=28., n=62, a=62.*g/mole);


  G4Element* Ni64enriched95 = new G4Element(name="Ni64enriched95",symbole="64Ni_95", ncomponents=5);
  Ni64enriched95->AddIsotope(Ni64, abundance = 94.8*perCent);
  Ni64enriched95->AddIsotope(Ni58, abundance = 2.67*perCent);
  Ni64enriched95->AddIsotope(Ni60, abundance = 1.75*perCent);
  Ni64enriched95->AddIsotope(Ni61, abundance = 0.11*perCent);
  Ni64enriched95->AddIsotope(Ni62, abundance = 0.67*perCent);


  fTarget_Material = new G4Material("FTarget_Material",fDensity_target,ncomponents=1);
  fTarget_Material->AddElement(Ni64enriched95, fractionmass=1.);
  */

  fTarget_Material =  nist->FindOrBuildMaterial("G4_Ni");

  //FOIL MATERIAL
  fFoil_Material =  nist->FindOrBuildMaterial("G4_Al");
 
  
  ////////////////////////////////////////////////////////////////////
  /////////////////////////////PART1//////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////
  //////////////////////////LAYER PART 1//////////////////////////////
  ////////////////////////////////////////////////////////////////////


  //Create the (external) layer around the beam vacuum tube

  G4double layer1_length_PART1         = 98.9*mm;
  G4double layer1_innerRadius_PART1    = 11.5*mm;
  G4double layer1_outerRadius_PART1    = 16.*mm;
  G4double layer1_hz_PART1             = 0.5*layer1_length_PART1;
 
  G4Tubs* solidLayer1_PART1 
    = new G4Tubs("Layer1_PART1",
		 layer1_innerRadius_PART1,
		 layer1_outerRadius_PART1,
		 layer1_hz_PART1,
		 startAngle,
		 spanningAngle);
 

  G4double layer2_length_PART1         = 124.6*mm;
  G4double layer2_innerRadius_PART1    = 7.5*mm;
  G4double layer2_outerRadius_PART1    = 11.5*mm;
  G4double layer2_hz_PART1             = 0.5*layer2_length_PART1;
 
  G4Tubs* solidLayer2_PART1 
    = new G4Tubs("Layer2_PART1",
		 layer2_innerRadius_PART1,
		 layer2_outerRadius_PART1,
		 layer2_hz_PART1,
		 startAngle,
		 spanningAngle);
  

 
  G4RotationMatrix rot_layer_PART1;
  G4double z_layer_translation_PART1 = layer2_length_PART1-layer1_length_PART1;
  G4ThreeVector placement_layer_PART1 = G4ThreeVector(0.*mm,0.*mm,0.5*z_layer_translation_PART1);
  G4Transform3D transform_layer_PART1(rot_layer_PART1,placement_layer_PART1);

  G4UnionSolid* solidLayer_PART1 = new G4UnionSolid("Layer_PART1", solidLayer2_PART1, solidLayer1_PART1, transform_layer_PART1);

  
  G4LogicalVolume* logicLayer_PART1
    = new G4LogicalVolume(solidLayer_PART1,
			  al,
			  "Layer_PART1");

  /*  G4VPhysicalVolume* physLayer1_PART1 = */
  new G4PVPlacement(0, 
		    G4ThreeVector(0.*mm,0.*mm,0.5*layer2_length_PART1),      
		    logicLayer_PART1,                                  
		    "Layer_PART1",                                    
		    fLogicWorld,                                       
		    false,                                            
		    0,                                                
		    checkOverlaps);                                    

  
  //////////////////////////////////////////////////////////////////
  ////////////////Create the beam vacuum tube///////////////////////
  //////////////////////////////////////////////////////////////////


  G4double tube_length_PART1 = layer2_length_PART1;
  
  G4double tube_innerRadius_PART1    = 0.*mm;
  G4double tube_outerRadius_PART1    = 7.5*mm;
  G4double tube_hz_PART1             = 0.5*tube_length_PART1;
  
  G4Tubs* solidTube_PART1 
    = new G4Tubs("Tube_PART1",
		 tube_innerRadius_PART1,
		 tube_outerRadius_PART1,
		 tube_hz_PART1,
		 startAngle,
		 spanningAngle);
  
  G4LogicalVolume* logicTube_PART1
    = new G4LogicalVolume(solidTube_PART1,
			  vacuum_beam,
			  "Tube_PART1");

  /*  G4VPhysicalVolume* physTube_PART1 = */
  new G4PVPlacement(0,                                               
		    G4ThreeVector(0.*mm,0.*mm,0.5*tube_length_PART1),
		    logicTube_PART1,                                 
		    "Tube_PART1",                                    
		    fLogicWorld,                               
		    false,                                           
		    0,                                              
		    checkOverlaps);                                  




  

  ////////////////////////////////////////////////////////////////////////////////
  /////////////////////////// PART 2 = FOIL PART /////////////////////////////////  
  ////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////
  //////////////////////////LAYER PART 2//////////////////////////////
  ////////////////////////////////////////////////////////////////////

  //Create the (external) layer around the foil part

  G4double layer_length_PART2          = 12.*mm;
  G4double layer_innerRadius_PART2     = 7.5*mm;
  G4double layer_outerRadius_PART2     = 15.*mm;
  G4double layer_hz_PART2              = 0.5*layer_length_PART2;
   
  G4Tubs* solidLayer_PART2 
    = new G4Tubs("Layer_PART2",
		 layer_innerRadius_PART2,
		 layer_outerRadius_PART2,
		 layer_hz_PART2,
		 startAngle,
		 spanningAngle);  

  
   G4LogicalVolume* logicLayer_PART2
    = new G4LogicalVolume(solidLayer_PART2,
			  al,
			  "Layer_PART2");

   G4double layer_z_position_PART2 = tube_length_PART1 + 0.5*layer_length_PART2;

   /*   G4VPhysicalVolume* physLayer1_PART2 = */ 
   new G4PVPlacement(0,                                                     
		     G4ThreeVector(0.*mm,0.*mm,layer_z_position_PART2),  
		     logicLayer_PART2,                                   
		     "Layer_PART2",                                     
		     fLogicWorld,                                         
		     false,                                              
		     0,                                                 
		     checkOverlaps);                                   



  ////////////////////////////////////////////////////////////////////
  //////////////////////////TUBE PART 2///////////////////////////////
  ////////////////////////////////////////////////////////////////////

  //Create the tube before the foil
     
  G4double tube_length_PART2 = layer_length_PART2;

  G4double tube_innerRadius_PART2     = 0.*mm;
  G4double tube_outerRadius_PART2     = 7.5*mm;
  G4double tube_hz_PART2              = 0.5*tube_length_PART2;
  
  G4Tubs* solidTube_PART2 
    = new G4Tubs("Tube_PART2",
		 tube_innerRadius_PART2,
		 tube_outerRadius_PART2,
		 tube_hz_PART2,
		 startAngle,
		 spanningAngle);
  
  G4LogicalVolume* logicTube_PART2
   = new G4LogicalVolume(solidTube_PART2,
			 vacuum_beam,
			 "Tube_PART2");

  
  /*  G4VPhysicalVolume* physTube_PART2 = */
  new G4PVPlacement(0,                                                 
		    G4ThreeVector(0.*mm,0.*mm, layer_z_position_PART2),
		    logicTube_PART2,                                   
		    "Tube_PART2",                                      
		    fLogicWorld,                                  
		    false,                                             
		    0,                                                
		    checkOverlaps);                                    
  

  /////////////////////////////////////////////////////////
  //////////////////////// GRID PART //////////////////////
  /////////////////////////////////////////////////////////


  //HEXAGONE
  G4double hexagone_length = layer_length_PART2 ;
  G4int numSide = 6;
  G4int numZPlanes = 2;
  const G4double zPlane[]={-0.5*hexagone_length,0.5*hexagone_length};
  const G4double rInner[]={3.*mm,3.*mm};
  const G4double rOuter[]={3.3*mm,3.3*mm};

  G4Polyhedra* solidHexagone_0
    = new G4Polyhedra("Grid",
		      0.*deg,
		      360.*deg,
		      numSide,
		      numZPlanes,
		      zPlane,
		      rInner,
		      rOuter);


  G4LogicalVolume* logicHexagone
    = new G4LogicalVolume(solidHexagone_0,
  			  al,
  			  "Grid");

  /*  G4VPhysicalVolume* physHexagone = */
  new G4PVPlacement(0,                                  
		   G4ThreeVector(0.*mm,0.*mm,0.*mm), 
		   logicHexagone,                    
		   "Grid",                           
		   logicTube_PART2,                  
		   false,                            
		   0,                                 
		   checkOverlaps);                    
 


 
 
  /////////////////////////////////////////////////////////////////////////
  //////////////////////////CREATE THE FOIL ///////////////////////////////
  /////////////////////////////////////////////////////////////////////////  

 
  fZ_foil_position = 0.5*fFoil_thickness + tube_length_PART1 + layer_length_PART2;
      

   
  G4double foil_innerRadius    = 0.*mm;
  G4double foil_outerRadius    = 16.*mm;
  
  fSolidFoil 
    = new G4Tubs("Foil",
		 foil_innerRadius,
		 foil_outerRadius,
		 0.5*fFoil_thickness,
		 startAngle,
		 spanningAngle);
  
  fFoilVolume = fSolidFoil->GetCubicVolume();
  
  fLogicFoil
    = new G4LogicalVolume(fSolidFoil,
			  fFoil_Material,
			  "Foil");
  
  
  
  fPhysFoil
    = new G4PVPlacement(0,                                    
			G4ThreeVector(0.*mm,0.*mm,fZ_foil_position),   
			fLogicFoil,                           
			"Foil",                              
			fLogicWorld,                          
			false,                               
			0,                                    
			checkOverlaps);                       
  




  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////// PART 3 = AFTER THE FOIL /////////////////////////////////  
  /////////////////////////////////////////////////////////////////////////////////////
  

  ////////////////////////////////////////////////////////////////////
  //////////////////////////LAYER PART 3//////////////////////////////
  ////////////////////////////////////////////////////////////////////

  //Create the (external) layer around the beam
  
  
  G4double layer_length_PART3 = 38.32*mm;

  G4double layer_innerRadius_PART3     = 7.5*mm;
  G4double layer_outerRadius_PART3     = 16.*mm;
  G4double layer_hz_PART3              = 0.5*layer_length_PART3;
  
  G4Tubs* solidLayer_PART3 
    = new G4Tubs("Layer_PART3",
		 layer_innerRadius_PART3,
		 layer_outerRadius_PART3,
		 layer_hz_PART3,
		 startAngle,
		 spanningAngle);
  
  G4LogicalVolume* logicLayer_PART3
    = new G4LogicalVolume(solidLayer_PART3,
			  al,
			  "Layer_PART3");

  fLayer_z_position_PART3 = tube_length_PART1 + layer_length_PART2 + fFoil_thickness + 0.5*layer_length_PART3;

  fPhysLayer_PART3
    = new G4PVPlacement(0,                                                  
  			G4ThreeVector(0.*mm,0.*mm,fLayer_z_position_PART3), 
  			logicLayer_PART3,                                   
  			"Layer_PART3",                                     
  			fLogicWorld,                                        
  			false,                                            
  			0,                                                  
  			checkOverlaps);                                     

  

  ////////////////////////////////////////////////////////////////////
  //////////////////////////TUBE PART 3///////////////////////////////
  ////////////////////////////////////////////////////////////////////

  //Create the tube after the foil, MATERIAL = HELIUM

  G4double tube_length_PART3 = layer_length_PART3;

  G4double tube_innerRadius_PART3     = 0.*mm;
  G4double tube_outerRadius_PART3     = 7.5*mm;
  G4double tube_hz_PART3              = 0.5*tube_length_PART3;
  
  G4Tubs* solidTube_PART3 
    = new G4Tubs("Tube_PART3",
		 tube_innerRadius_PART3,
		 tube_outerRadius_PART3,
		 tube_hz_PART3,
		 startAngle,
		 spanningAngle);
  
  G4LogicalVolume* logicTube_PART3
    = new G4LogicalVolume(solidTube_PART3,
			  helium,
			  "Tube_PART3");

 
  fPhysTube_PART3
    = new G4PVPlacement(0,                                                  
			G4ThreeVector(0.*mm,0.*mm,fLayer_z_position_PART3), 
			logicTube_PART3,                                   
			"Tube_PART3",                                      
			fLogicWorld,                                        
			false,                                             
			0,                                                  
			checkOverlaps);                                    


  ////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////// PART 4 = BEFORE THE TARGET/////////////////////////////////  
  ////////////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////
  //////////////////////////LAYER PART4///////////////////////////////
  ////////////////////////////////////////////////////////////////////  

 
  G4double layer1_length_PART4    = 11.5*mm; //Length of the gold part

  G4double layer2_length_PART4    = 4.6*mm;
  G4double layer2_Rmin1_PART4     = 8.5*mm;
  G4double layer2_Rmax1_PART4     = 9.*mm;
  G4double layer2_Rmin2_PART4     = 8.5*mm;
  G4double layer2_Rmax2_PART4     = 14.*mm;
  G4double layer2_hz_PART4        = 0.5*layer2_length_PART4;
  
  G4Cons* solidLayer2_PART4 
    = new G4Cons("Layer2_PART4",
		 layer2_Rmin1_PART4,
		 layer2_Rmax1_PART4,
		 layer2_Rmin2_PART4,
		 layer2_Rmax2_PART4,
		 layer2_hz_PART4,
		 startAngle,
		 spanningAngle);
  
  //Create the (external) aluminium layer around the beam before the target

  G4double layer3_length_PART4 = layer1_length_PART4-layer2_length_PART4;
  G4double layer3_innerRadius_PART4     = 8.5*mm;
  G4double layer3_outerRadius_PART4     = 14.*mm;
  G4double layer3_hz_PART4              = 0.5*layer3_length_PART4;
  
  G4Tubs* solidLayer3_PART4 
    = new G4Tubs("Layer3_PART4",
		 layer3_innerRadius_PART4,
		 layer3_outerRadius_PART4,
		 layer3_hz_PART4,
		 startAngle,
		 spanningAngle);
  
  

  G4RotationMatrix rot_layer_PART4;
  G4double z_layer_translation_PART4 = 0.5*(layer2_length_PART4+layer3_length_PART4);
  G4ThreeVector placement_layer_PART4 = G4ThreeVector(0.*mm,0.*mm,z_layer_translation_PART4);
  G4Transform3D transform_layer_PART4(rot_layer_PART4,placement_layer_PART4);

  G4UnionSolid* solidLayer_PART4 = new G4UnionSolid("Layer_PART1", solidLayer2_PART4, solidLayer3_PART4, transform_layer_PART4);



  G4LogicalVolume* logicLayer_PART4
    = new G4LogicalVolume(solidLayer_PART4,
			  al,
			  "Layer_PART4");

  fLayer_z_position_PART4 = tube_length_PART1 + layer_length_PART2 + fFoil_thickness + layer_length_PART3 +0.5*layer2_length_PART4;

  fPhysLayer_PART4
    = new G4PVPlacement(0,                                                  
			G4ThreeVector(0.*mm,0.*mm,fLayer_z_position_PART4),  
			logicLayer_PART4,                                  
			"Layer_PART4",                                    
			fLogicWorld,                                        
			false,                                              
			0,                                                 
			checkOverlaps);                                   





   //Create the (internal) gold layer around the beam

  G4double layer1_innerRadius_PART4     = 8.*mm;
  G4double layer1_outerRadius_PART4     = 8.5*mm;
  G4double layer1_hz_PART4              = 0.5*layer1_length_PART4;
  
  G4Tubs* solidLayer1_PART4 
    = new G4Tubs("Layer1_PART4",
		 layer1_innerRadius_PART4,
		 layer1_outerRadius_PART4,
		 layer1_hz_PART4,
		 startAngle,
		 spanningAngle);

 
  G4LogicalVolume* logicLayer1_PART4
    = new G4LogicalVolume(solidLayer1_PART4,
			  pt,
			  "Layer1_PART4");
  
  fLayer1_z_position_PART4 = tube_length_PART1 + layer_length_PART2 + fFoil_thickness + layer_length_PART3 +0.5*layer1_length_PART4;

  
  fPhysLayer1_PART4
    = new G4PVPlacement(0,                                                 
			G4ThreeVector(0.*mm,0.*mm,fLayer1_z_position_PART4),
			logicLayer1_PART4,                                 
			"Layer1_PART4",                                   
			fLogicWorld,                                  
			false,                                     
			0,                              
			checkOverlaps);              

  ////////////////////////////////////////////////////////////////////
  //////////////////////////TUBE PART 4///////////////////////////////
  ////////////////////////////////////////////////////////////////////


  //Create the tube before the target

 
  fTube_length_PART4 = layer1_length_PART4;

  G4double tube_innerRadius_PART4     = 0.*mm;
  fTube_outerRadius_PART4     = 8.*mm;
  G4double tube_hz_PART4              = 0.5*fTube_length_PART4;
 
  G4Tubs* solidTube_PART4 
    = new G4Tubs("Tube_PART4",
		 tube_innerRadius_PART4,
		 fTube_outerRadius_PART4,
		 tube_hz_PART4,
		 startAngle,
		 spanningAngle);
  
  G4LogicalVolume* logicTube_PART4
    = new G4LogicalVolume(solidTube_PART4,
			  helium,
			  "Tube_PART4");


  fPhysTube_PART4
    = new G4PVPlacement(0,                                                  
			G4ThreeVector(0.*mm,0.*mm,fLayer1_z_position_PART4), 
			logicTube_PART4,                                    
			"Tube_PART4",                                       
			fLogicWorld,                                         
			false,                                              
			0,                                                  
			checkOverlaps);                                    




  ////////////////////////////////////////////////////////////////////
  //////////////////////////// TARGET ////////////////////////////////
  ////////////////////////////////////////////////////////////////////


  

  //Design the target inside the tube

  G4double target_innerRadius     = 0.*mm;
  G4double target_outerRadius     = 0.5*fTarget_diameter;
  G4double target_hz              = 0.5*fTarget_thickness;
 
  fSolidTarget 
    = new G4Tubs("Target",
		 target_innerRadius,
		 target_outerRadius,
		 target_hz,
		 startAngle,
		 spanningAngle);

  fTargetVolume = fSolidTarget->GetCubicVolume();
  
  fLogicTarget
    = new G4LogicalVolume(fSolidTarget,
			  fTarget_Material,
			  "Target");

  fTarget_z_position = 0.5*fTube_length_PART4-0.5*fTarget_thickness;  //inside the tube!


  fPhysTarget
    = new G4PVPlacement(0,                                                  
			G4ThreeVector(0.*mm,0.*mm, fTarget_z_position),      
			fLogicTarget,                                       
			"Target",                                          
			logicTube_PART4,                                 
			false,                                           
			0,                                               
			checkOverlaps);                                   

 

  ///////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////// PART 5 = AFTER THE TARGET /////////////////////////////////  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////
  //////////////////////////LAYER PART 5//////////////////////////////
  ////////////////////////////////////////////////////////////////////

  //Create the (internal) platinium layer

  G4double layer1_length_PART5 = 0.5*mm;

  G4double layer1_innerRadius_PART5     = 0.*mm;
  G4double layer1_outerRadius_PART5     = 8.5*mm;
  G4double layer1_hz_PART5              = 0.5*layer1_length_PART5;
  
  G4Tubs* solidLayer1_PART5 
    = new G4Tubs("Layer1_PART5",
		 layer1_innerRadius_PART5,
		 layer1_outerRadius_PART5,
		 layer1_hz_PART5,
		 startAngle,
		 spanningAngle);
  
  G4LogicalVolume* logicLayer1_PART5
    = new G4LogicalVolume(solidLayer1_PART5,
			  pt,
			  "Layer1_PART5");

  fLayer1_z_position_PART5 = tube_length_PART1 + layer_length_PART2 + fFoil_thickness + layer_length_PART3 +layer1_length_PART4 + 0.5*layer1_length_PART5;

  fPhysLayer1_PART5
    = new G4PVPlacement(0,                                                  
			G4ThreeVector(0.*mm,0.*mm,fLayer1_z_position_PART5), 
			logicLayer1_PART5,                                 
			"Layer1_PART5",                                   
			fLogicWorld,                                     
			false,                                      
			0,                                              
			checkOverlaps);                                  
 


  //Create the (external) aluminium layer after the target

  G4double layer2_length_PART5 = 0.5*mm;

  G4double layer2_innerRadius_PART5     = 8.5*mm;
  G4double layer2_outerRadius_PART5     = 14.*mm;
  G4double layer2_hz_PART5              = 0.5*layer2_length_PART5;
  
  G4Tubs* solidLayer2_PART5 
    = new G4Tubs("Layer2_PART5",
		 layer2_innerRadius_PART5,
		 layer2_outerRadius_PART5,
		 layer2_hz_PART5,
		 startAngle,
		 spanningAngle);
  
  G4LogicalVolume* logicLayer2_PART5
    = new G4LogicalVolume(solidLayer2_PART5,
			  al,
			  "Layer2_PART5");

  fLayer2_z_position_PART5 = tube_length_PART1 + layer_length_PART2 + fFoil_thickness + layer_length_PART3 +layer2_length_PART4+layer3_length_PART4 + 0.5*layer2_length_PART5;

  fPhysLayer2_PART5
    = new G4PVPlacement(0,                                                   
			G4ThreeVector(0.*mm,0.*mm,fLayer2_z_position_PART5),
			logicLayer2_PART5,                                 
			"Layer2_PART5",                                     
			fLogicWorld,                                         
			false,                                              
			0,                                                 
			checkOverlaps);                                    

 

  //Create the end of the beam target machine

  G4double layer3_length_PART5 = 3.*mm;

  G4double layer3_innerRadius_PART5     = 0.*mm;
  G4double layer3_outerRadius_PART5     = 14.*mm;
  G4double layer3_hz_PART5              = 0.5*layer3_length_PART5;
  
  G4Tubs* solidLayer3_PART5 
    = new G4Tubs("Layer3_PART5",
		 layer3_innerRadius_PART5,
		 layer3_outerRadius_PART5,
		 layer3_hz_PART5,
		 startAngle,
		 spanningAngle);
  
  G4LogicalVolume* logicLayer3_PART5
    = new G4LogicalVolume(solidLayer3_PART5,
			  al,
			  "Layer3_PART5");

  fLayer3_z_position_PART5 = tube_length_PART1 + layer_length_PART2 + fFoil_thickness + layer_length_PART3 +layer2_length_PART4+layer3_length_PART4 + layer2_length_PART5 + 0.5*layer3_length_PART5;

  fPhysLayer3_PART5
    = new G4PVPlacement(0,                                                   
			G4ThreeVector(0.*mm,0.*mm,fLayer3_z_position_PART5), 
			logicLayer3_PART5,                                  
			"Layer3_PART5",                                     
			fLogicWorld,                                        
			false,                                             
			0,                                                  
			checkOverlaps);                                  




  //////////////////////////////////////////////
  /////////// Set sensitive region /////////////
  //////////////////////////////////////////////


  if(!fRegionTarget)
    {
      fRegionTarget = new G4Region("Target");
      fLogicTarget -> SetRegion(fRegionTarget);
      fRegionTarget -> AddRootLogicalVolume(fLogicTarget);
    }

  if(!fRegionFoil&&fPhysFoil!=nullptr)
    {
      fRegionFoil = new G4Region("Foil");
      fLogicFoil -> SetRegion(fRegionFoil);
      fRegionFoil -> AddRootLogicalVolume(fLogicFoil);
    }
  
  
 

  //Always return physical world//

  return physWorld;


}

void STCyclotronDetectorConstruction::ConstructSDandField()
{

  if(fLogicTarget != nullptr){
    STCyclotronSensitiveTarget* TargetSD = new STCyclotronSensitiveTarget("TargetSD", this);
    SetSensitiveDetector("Target",TargetSD);
  }

  
  if(fLogicFoil != nullptr){
    STCyclotronSensitiveFoil* FoilSD = new STCyclotronSensitiveFoil("FoilSD", this);
    SetSensitiveDetector("Foil",FoilSD);
  }
}

void STCyclotronDetectorConstruction::SetTargetDiameter(G4double targetDiameter)
{

  if(fTarget_diameter!=targetDiameter){
    if(targetDiameter/2.>fTube_outerRadius_PART4){ G4cout << "Error : the diameter is bigger than the tube" << G4endl;}
    else{
      fTarget_diameter = targetDiameter;
      if(fSolidTarget) fSolidTarget->SetOuterRadius(0.5*fTarget_diameter);
      G4cout << "The new diameter of the target is " << fTarget_diameter << " mm." << G4endl;
    }
  }

  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
  
}

//SET MATERIAL METHODS//

void STCyclotronDetectorConstruction::SetTargetIsotopeName(G4String name){
  fIsotopeName.push_back(name);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetTargetIsotopeZ(G4double z){
  fIsotopeZ.push_back(z);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetTargetIsotopeN(G4int n){
  fIsotopeN.push_back(n);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetTargetIsotopeA(G4double a){
  fIsotopeA.push_back(a);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetTargetElementName(G4String name){
  fElementName.push_back(name);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetTargetElementSymbole(G4String symbole){
  fElementSymbole.push_back(symbole);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetTargetElementNComponents(G4int n){
  fElementNComponents.push_back(n);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetTargetElementAbundance(G4double abundance){
  fElementAbundance.push_back(abundance);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetTargetMaterialDensity(G4double materialDensity){
  fDensity_target = materialDensity;
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetTargetMaterialNComponents(G4int n){
  fTarget_NComponents = n;
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetTargetMaterialFractionMass(G4double fractionMass){
  fMaterialFractionMass.push_back(fractionMass);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetTargetNaturalElement(G4String name){
  fNaturalElementName.push_back(name);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetTargetNaturalMaterialFractionMass(G4double fractionMass){
  fNaturalMaterialFractionMass.push_back(fractionMass);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

//UPDATE THE MATERIAL TO APPLY THE MODIFICATIONS

G4bool STCyclotronDetectorConstruction::UpdateMaterial(){
  
  G4int nElements = fTarget_NComponents;

  G4double density = fDensity_target*g/cm3;
  G4Material* material = new G4Material("FTarget_Material", density, nElements);
  
  G4double checkFractionMass = 0;

 
  //CHECK THAT NUMBER OF ELEMENTS IN THE MATERIAL IS EQUAL TO THE NUMBER OF ELEMENTS DEFINED USING ISOTOPES PLUS THE NUMBER OF ELEMENT DEFINED USING NIST
 
  G4int counterElement = 0;
  
  //std::vector<G4String>::size_type sz = fElementName.size();
  std::ptrdiff_t const sizeElementName = fElementName.size();
  std::ptrdiff_t const sizeNaturalElementName = fNaturalElementName.size();


  if(sizeElementName!=0){

    //ELEMENT LOOP
    if(nElements != sizeElementName + sizeNaturalElementName){
      G4cout << "Error : the number of elements in the target was set up at " << nElements << " but you defined only " << fElementName.size()+fNaturalElementName.size() << "elements." << G4endl ;
      return false;
    }
    
    for(int i = 0; i<sizeElementName ; i++){
      
      checkFractionMass = checkFractionMass + fMaterialFractionMass[i];
      G4Element* elementi = new G4Element(fElementName[i],fElementSymbole[i],fElementNComponents[i]);
      
      G4double checkAbundance = 0.;
    
      //ISOTOPE LOOP FOR THE ELEMENT
      std::ptrdiff_t const sizeIsotopeName = fIsotopeName.size();

      if(fElementNComponents[i] != sizeIsotopeName){
	G4cout << "Error : the number of isotopes defined in the target's element" << fElementName[i] << "was set up at " << fElementNComponents[i] << " but you defined only " << fIsotopeName.size() << "elements." << G4endl;
	return false;
      }
      for(G4int j=counterElement;j<counterElement+fElementNComponents[i];++j){
	G4double A = fIsotopeA[j]*g/mole;
	G4Isotope* isotopeij = new G4Isotope(fIsotopeName[j], G4int(fIsotopeZ[j]), fIsotopeN[j], A);
	checkAbundance = checkAbundance + fElementAbundance[j];
	elementi->AddIsotope(isotopeij,fElementAbundance[j]);
      }
      if((1.-checkAbundance)>1E-5){
	G4cout << "Error : the total abundance of isotopes in your target's element " << fElementName[i] << " is equal to " << checkAbundance << ". It must be equal to 1." << G4endl;
	return false;
      }

      counterElement = counterElement + fElementNComponents[i];
      material->AddElement(elementi, fMaterialFractionMass[i]);
      
    }
  }


  if(sizeNaturalElementName!=0){
    for(int i=0;i<sizeNaturalElementName;i++){
      checkFractionMass = checkFractionMass + fNaturalMaterialFractionMass[i];
      G4NistManager* man = G4NistManager::Instance();
      G4Element* element = man->FindOrBuildElement(fNaturalElementName[i]);
      material->AddElement(element, fNaturalMaterialFractionMass[i]);
    }
  }
  
  if((1.-checkFractionMass)>1E-5){
    G4cout << "Error : the sum of the fraction mass of each element in the target is equal to " << checkFractionMass << ". It must be equal to 1." << G4endl;
    return false;
  }

  fTarget_Material = material;
  G4cout << "You succesfully changed the material of the target." << G4endl;
  G4cout << "Here is the new material : " << G4endl;
  G4cout << "Number of elements : " << material->GetNumberOfElements() << G4endl;

  std::ptrdiff_t const sizeMaterial = material->GetNumberOfElements();
 

  for(int i = 0; i<sizeMaterial;i++){
    const G4Element* element = material->GetElement(i);
    G4cout << element->GetName() << " having the following isotopes : " << G4endl;

     std::ptrdiff_t const sizeIsotope = element->GetNumberOfIsotopes();
 

    for(int j=0; j<sizeIsotope; j++){
      const G4Isotope* isotope = element->GetIsotope(j);
      G4cout << " - " << (element->GetRelativeAbundanceVector())[j]*100 << "% of " << isotope->GetName() << " Z = " << isotope->GetZ() << " N = " << isotope->GetN() << "." << G4endl;
    }
  }
  fLogicTarget->SetMaterial(fTarget_Material);

  //
  fIsotopeName.clear();
  fIsotopeZ.clear();
  fIsotopeN.clear();
  fIsotopeA.clear();
  fElementName.clear();
  fElementSymbole.clear();
  fElementNComponents.clear();
  fElementAbundance.clear();
  fNaturalElementName.clear();
  fNaturalMaterialFractionMass.clear();
  fMaterialFractionMass.clear();
  //

  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
  
  
  return true;

  

}

//SET TARGET MATERIAL WITH PHYSICS NIST LIST//

void STCyclotronDetectorConstruction::SetTargetMaterial(G4String materialName){

  G4NistManager* nistManager = G4NistManager::Instance();
  G4Material* targetMaterial = nistManager->FindOrBuildMaterial(materialName);
  if(fTarget_Material!=targetMaterial){
    
    if(targetMaterial){
      
    fTarget_Material = targetMaterial;
    fLogicTarget->SetMaterial(fTarget_Material);
    
    G4cout << "The new material of the target is : " << fTarget_Material << G4endl;

    }
    else{ G4cout << "This material wasn't found in the NIST list." << G4endl; }
    
  }

  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;

}


void STCyclotronDetectorConstruction::SetFoilIsotopeName(G4String name){
  fIsotopeNameFoil.push_back(name);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetFoilIsotopeZ(G4double z){
  fIsotopeZFoil.push_back(z);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetFoilIsotopeN(G4int n){
  fIsotopeNFoil.push_back(n);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetFoilIsotopeA(G4double a){
  fIsotopeAFoil.push_back(a);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetFoilElementName(G4String name){
  fElementNameFoil.push_back(name);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetFoilElementSymbole(G4String symbole){
  fElementSymboleFoil.push_back(symbole);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetFoilElementNComponents(G4int n){
  fElementNComponentsFoil.push_back(n);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetFoilElementAbundance(G4double abundance){
  fElementAbundanceFoil.push_back(abundance);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetFoilMaterialDensity(G4double materialDensity){
  fDensity_foil = materialDensity;
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetFoilMaterialNComponents(G4int n){
  fFoil_NComponents = n;
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetFoilMaterialFractionMass(G4double fractionMass){
  fMaterialFractionMassFoil.push_back(fractionMass);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetFoilNaturalElement(G4String name){
  fNaturalElementNameFoil.push_back(name);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

void STCyclotronDetectorConstruction::SetFoilNaturalMaterialFractionMass(G4double fractionMass){
  fNaturalMaterialFractionMassFoil.push_back(fractionMass);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
}

G4bool STCyclotronDetectorConstruction::UpdateFoilMaterial(){
  
  G4int nElements = fFoil_NComponents;
  G4double density = fDensity_foil*g/cm3;
  
  G4Material* material = new G4Material("FFoil_Material", density,nElements);
  
  G4double checkFractionMass = 0;

 
  //CHECK THAT NUMBER OF ELEMENTS IN THE MATERIAL IS EQUAL TO THE NUMBER OF ELEMENTS DEFINED USING ISOTOPES PLUS THE NUMBER OF ELEMENT DEFINED USING NIST
 
  G4int counterElement = 0;
 
 std::ptrdiff_t const sizeElementName = fElementNameFoil.size();
  std::ptrdiff_t const sizeNaturalElementName = fNaturalElementNameFoil.size();
  

 
  if(sizeElementName!=0){
    //ELEMENT LOOP
    if(nElements != sizeElementName + sizeNaturalElementName){
      G4cout << "Error : the number of elements was set up at " << nElements << " but you defined only " << fElementNameFoil.size()+fNaturalElementNameFoil.size() << "elements." << G4endl ;
      return false;
    }
    for(int i = 0; i<sizeElementName ; i++){
      
      checkFractionMass = checkFractionMass + fMaterialFractionMassFoil[i];
      G4Element* elementi = new G4Element(fElementNameFoil[i],fElementSymboleFoil[i],fElementNComponentsFoil[i]);
      
      G4double checkAbundance = 0.;
      //ISOTOPE LOOP FOR THE ELEMENT

      std::ptrdiff_t const sizeIsotopeName = fIsotopeNameFoil.size();

      if(fElementNComponentsFoil[i] != sizeIsotopeName){
	G4cout << "Error : the number of isotopes in element" << fElementNameFoil[i] << " of the foil was set up at " << fElementNComponentsFoil[i] << " but you defined only " << fIsotopeNameFoil.size() << "elements." << G4endl;
	return false;
      }
      for(G4int j=counterElement;j<counterElement+fElementNComponentsFoil[i];++j){
	G4double A = fIsotopeAFoil[j]*g/cm3;
	G4Isotope* isotopeij = new G4Isotope(fIsotopeNameFoil[j], G4int(fIsotopeZFoil[j]), fIsotopeNFoil[j],A);
	checkAbundance = checkAbundance + fElementAbundanceFoil[j];
	elementi->AddIsotope(isotopeij,fElementAbundanceFoil[j]);
      }
      if((1.-checkAbundance)>1E-5){
	G4cout << "Error : the total abundance of isotopes in your foil's element " << fElementNameFoil[i] << " is equal to " << checkAbundance << ". It must be equal to 1." << G4endl;
	return false;
      }

      counterElement = counterElement + fElementNComponentsFoil[i];
      material->AddElement(elementi, fMaterialFractionMassFoil[i]);
      
    }
  }

  if(sizeNaturalElementName!=0){
    for(int i=0;i<sizeNaturalElementName;i++){
      checkFractionMass = checkFractionMass + fNaturalMaterialFractionMassFoil[i];
      G4NistManager* man = G4NistManager::Instance();
      G4Element* element = man->FindOrBuildElement(fNaturalElementNameFoil[i]);
      material->AddElement(element, fNaturalMaterialFractionMassFoil[i]);
    }
  }
  
  if((1.-checkFractionMass)>1E-5){
    G4cout << "Error : the sum of the fraction mass of each element in the foil is equal to " << checkFractionMass << ". It must be equal to 1." << G4endl;
    return false;
  }

  fFoil_Material = material;
  G4cout << "You succesfully changed the material of the foil." << G4endl;
  G4cout << "Here is the new material : " << G4endl;
  G4cout << "Number of elements : " << material->GetNumberOfElements() << G4endl;
  std::ptrdiff_t const sizeMaterial = material->GetNumberOfElements();


  for(int i = 0; i<sizeMaterial;i++){
    const G4Element* element = material->GetElement(i);
    G4cout << element->GetName() << " having the following isotopes : " << G4endl;
    
    std::ptrdiff_t const sizeIsotope = element->GetNumberOfIsotopes();
    
    for(int j=0; j<sizeIsotope; j++){
      const G4Isotope* isotope = element->GetIsotope(j);
      G4cout << " - " << (element->GetRelativeAbundanceVector())[j]*100 << "% of " << isotope->GetName() << " Z = " << isotope->GetZ() << " N = " << isotope->GetN() << "." << G4endl;
    }
  }
  fLogicFoil->SetMaterial(fFoil_Material);

  //
  fIsotopeNameFoil.clear();
  fIsotopeZFoil.clear();
  fIsotopeNFoil.clear();
  fIsotopeAFoil.clear();
  fElementNameFoil.clear();
  fElementSymboleFoil.clear();
  fElementNComponentsFoil.clear();
  fElementAbundanceFoil.clear();
  fNaturalElementNameFoil.clear();
  fNaturalMaterialFractionMassFoil.clear();
  fMaterialFractionMassFoil.clear();
  //

  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;

  return true;
}


void STCyclotronDetectorConstruction::SetFoilMaterial(G4String materialName){

  G4NistManager* nistManager = G4NistManager::Instance();
  G4Material* foilMaterial = nistManager->FindOrBuildMaterial(materialName);
  if(fFoil_Material!=foilMaterial){
    
    if(foilMaterial){
      
    fFoil_Material = foilMaterial;
    fLogicFoil->SetMaterial(fFoil_Material);
    
    G4cout << "The new material of the foil is : " << fFoil_Material << G4endl;

    }
    else{ G4cout << "This material wasn't found in the NIST list." << G4endl; }
    
  }
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;

}

//SET PARAMETERS FOR THE TARGET

void STCyclotronDetectorConstruction::SetTargetThickness(G4double targetThickness)
{

   if(fTarget_thickness!=targetThickness){
    if(targetThickness > fTube_length_PART4){ G4cout << "Error : the target thickness is longer than the tube length" << G4endl; }
    else {
      //Modify the position of the target
      fTarget_z_position = fTarget_z_position + 0.5*fTarget_thickness;
      fTarget_thickness = targetThickness;
      fTarget_z_position = fTarget_z_position - 0.5*fTarget_thickness;
      fPhysTarget->SetTranslation(G4ThreeVector(0.*mm,0.*mm,fTarget_z_position*mm));
      //Modify the thickness of the target
      fSolidTarget->SetZHalfLength(0.5*targetThickness);
      G4cout << "The new thickness of the target is " << fTarget_thickness << " mm." << G4endl;
      G4cout << "Position 1 = " << GetTargetPosition1() << " -- Position 2 = " << GetTargetPosition2() << G4endl;

    }
  }

   G4RunManager::GetRunManager() -> GeometryHasBeenModified();
   G4cout << "... Geometry is notified .... " << G4endl;
}


void STCyclotronDetectorConstruction::SetFoilThickness(G4double foilThickness)
{
  
  if(fFoil_thickness != foilThickness){
    
    
    //Change de position of the detector parts after the foil
    //Change the position to set it so there is no foil
    fLayer_z_position_PART3 = fLayer_z_position_PART3 - fFoil_thickness;
    fLayer_z_position_PART4 = fLayer_z_position_PART4 - fFoil_thickness;
    fLayer1_z_position_PART4 = fLayer1_z_position_PART4 - fFoil_thickness;
    fLayer1_z_position_PART5 = fLayer1_z_position_PART5 - fFoil_thickness;
    fLayer2_z_position_PART5 = fLayer2_z_position_PART5 - fFoil_thickness;
    fLayer3_z_position_PART5 = fLayer3_z_position_PART5 - fFoil_thickness;
    fZ_foil_position = fZ_foil_position - 0.5*fFoil_thickness;
    
    fFoil_thickness = foilThickness;
    
    //Change the position using the new foil thickness
    fLayer_z_position_PART3 = fLayer_z_position_PART3 + fFoil_thickness;
    fLayer_z_position_PART4 = fLayer_z_position_PART4 + fFoil_thickness;
    fLayer1_z_position_PART4 = fLayer1_z_position_PART4 + fFoil_thickness;
    fLayer1_z_position_PART5 = fLayer1_z_position_PART5 + fFoil_thickness;
    fLayer2_z_position_PART5 = fLayer2_z_position_PART5 + fFoil_thickness;
    fLayer3_z_position_PART5 = fLayer3_z_position_PART5 + fFoil_thickness;
    fZ_foil_position = fZ_foil_position + 0.5*fFoil_thickness;
    
    fPhysLayer_PART3->SetTranslation(G4ThreeVector(0.*mm,0.*mm,fLayer_z_position_PART3*mm));
    fPhysTube_PART3->SetTranslation(G4ThreeVector(0.*mm,0.*mm,fLayer_z_position_PART3*mm));
    fPhysLayer_PART4->SetTranslation(G4ThreeVector(0.*mm,0.*mm,fLayer_z_position_PART4*mm));
    fPhysLayer1_PART4->SetTranslation(G4ThreeVector(0.*mm,0.*mm,fLayer1_z_position_PART4*mm));
    fPhysTube_PART4->SetTranslation(G4ThreeVector(0.*mm,0.*mm,fLayer1_z_position_PART4*mm));
    fPhysLayer1_PART5->SetTranslation(G4ThreeVector(0.*mm,0.*mm,fLayer1_z_position_PART5*mm));
    fPhysLayer2_PART5->SetTranslation(G4ThreeVector(0.*mm,0.*mm,fLayer2_z_position_PART5*mm));
    fPhysLayer3_PART5->SetTranslation(G4ThreeVector(0.*mm,0.*mm,fLayer3_z_position_PART5*mm));
    
    
    if(fPhysFoil) fPhysFoil->SetTranslation(G4ThreeVector(0.*mm,0.*mm,fZ_foil_position*mm));
    //Change the thickness
    if(fSolidFoil) fSolidFoil->SetZHalfLength(0.5*foilThickness*mm);
    G4cout << "The new thickness of the foil is " << fFoil_thickness << " mm." << G4endl;
        
  }

  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... " << G4endl;
  
}


