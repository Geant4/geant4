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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// If you use this example, please cite the following publication:
// Rad. Prot. Dos. 133 (2009) 2-11

#include "DetectorConstruction.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreadLocal EMField* DetectorConstruction::fField = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::DetectorConstruction()
  
  :fDefaultMaterial(NULL),fCollimatorMaterial(NULL),fBoiteMaterial(NULL),
   fCathodeMaterial(NULL),fVerreMaterial(NULL),fVerre2Material(NULL),   
   fKgmMaterial(NULL),fBoite2Material(NULL),fBoite3Material(NULL),
   fNucleusMaterial1(NULL),fCytoplasmMaterial1(NULL),
   fNucleusMaterial2(NULL),fCytoplasmMaterial2(NULL),
   fNucleusMaterial3(NULL),fCytoplasmMaterial3(NULL),
   fPhysiWorld(NULL),fLogicWorld(NULL),fSolidWorld(NULL),
   fPhysiVol(NULL),fLogicVol(NULL),fSolidVol(NULL),
   fPhysiBoite(NULL),fLogicBoite(NULL),fSolidBoite(NULL),
   fPhysiYoke1(NULL),fLogicYoke1(NULL),fSolidYoke1(NULL),
   fPhysi1Gap(NULL),fLogic1Gap(NULL),fSolid1Gap(NULL),
   fPhysi2Gap(NULL),fLogic2Gap(NULL),fSolid2Gap(NULL), 
   fPhysi3Gap(NULL),fLogic3Gap(NULL),fSolid3Gap(NULL),
   fPhysiYoke2(NULL),fLogicYoke2(NULL),fSolidYoke2(NULL),
   fPhysi4Gap(NULL),fLogic4Gap(NULL),fSolid4Gap(NULL),
   fPhysi5Gap(NULL),fLogic5Gap(NULL),fSolid5Gap(NULL), 
   fPhysiBoiteIso(NULL),fLogicBoiteIso(NULL),fSolidBoiteIso(NULL),
   fPhysiCathode(NULL),fLogicCathode(NULL),fSolidCathode(NULL), 
   fPhysiIso(NULL),fLogicIso(NULL),fSolidIso(NULL),
   fPhysiVerre(NULL),fLogicVerre(NULL),fSolidVerre(NULL),
   fPhysiBoite2(NULL),fLogicBoite2(NULL),fSolidBoite2(NULL),
   fPhysiBoite3(NULL),fLogicBoite3(NULL),fSolidBoite3(NULL),
   fPhysiKgm(NULL),fLogicKgm(NULL),fSolidKgm(NULL),
   fPhysiVerre2(NULL),fLogicVerre2(NULL),fSolidVerre2(NULL),
   fPhysiPhantom(NULL),fLogicPhantom(NULL),fSolidPhantom(NULL)
  
{
  fWorldSizeXY=fWorldSizeZ=0;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::~DetectorConstruction()
{
  delete fDefaultMaterial;
  delete fCollimatorMaterial;
  delete fBoiteMaterial;
  delete fCathodeMaterial;
  delete fVerreMaterial;
  delete fVerre2Material;
  delete fKgmMaterial;
  delete fBoite2Material;
  delete fBoite3Material;
  delete fNucleusMaterial1;
  delete fCytoplasmMaterial1;
  delete fNucleusMaterial2;
  delete fCytoplasmMaterial2;
  delete fNucleusMaterial3;
  delete fCytoplasmMaterial3;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* DetectorConstruction::Construct()
  
{
  DefineMaterials();
  return ConstructLine();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::DefineMaterials()
{ 

  G4String name, symbol;             
  G4double density;            
  
  G4int ncomponents, natoms,nel;
  G4double z, a;
  G4double fractionmass;
  G4double temperature, pressure;
  
  // Define Elements 
  
  G4Element*   H  = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);
  G4Element*   N  = new G4Element ("Nitrogen", "N", 7., 14.01*g/mole);
  G4Element*   O  = new G4Element ("Oxygen"  , "O", 8. , 16.00*g/mole);
  G4Element*   Ar = new G4Element ("Argon" , "Ar", 18., 39.948*g/mole );
  G4Element*    C = new G4Element ("Carbon","C", 6., 12.011*g/mole);
  G4Element *  Si = new G4Element ("Silicon","Si",14., 28.0855*g/mole);
  G4Element *  Cu = new G4Element ("Cuivre","Cu",29., 63.546*g/mole);
  G4Element *  Zn = new G4Element ("Zinc","Zn",30.,65.409*g/mole);
  G4Element *  P  = new G4Element ("Phosphorus","P",15.,30.973761*g/mole);
 
  // Vacuum standard definition...
  
  density = universe_mean_density;
  G4Material* vacuum = new G4Material(name="Vacuum", z=1., a=1.01*g/mole,
				      density);  
  // Water 
  
  density = 1.000*g/cm3;
  G4Material* H2O = new G4Material(name="H2O"  , density, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
 
  // Air
  
  density = 1.290*mg/cm3;
  pressure = 1*atmosphere;
  temperature = 293.16*kelvin;
  G4Material* Air = new G4Material(name="Air"  , density, ncomponents=2, kStateGas, temperature, pressure);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);
  
  // Low Pressure air
  
  density = (5e-6/1013.)*1.290*mg/cm3; // 5e-6 mbar is the usual beam pipe air pressure
  pressure = 1*atmosphere;
  temperature = 293.16*kelvin;
  G4Material* LPAir = new G4Material(name="LPAir"  , density, ncomponents=3, kStateGas, temperature, pressure);
  LPAir->AddElement(N, fractionmass=0.715);
  LPAir->AddElement(O, fractionmass=0.25);
  LPAir->AddElement(Ar, fractionmass=0.035);
  
  // Platinum
  
  a = 195.09*g/mole;
  density = 21.4*g/cm3;
  G4Material* Pt = new G4Material(name="Pl", z=78., a, density);
 
  // Butane @ 10 mbar
  
  density = 2.552e-2*mg/cm3;
  pressure = 0.01*bar;
  temperature = 293.16*kelvin;
  G4Material* Butane = new G4Material(name = "Butane", density, nel = 2, kStateGas, temperature, pressure);
  Butane->AddElement (C, natoms=4);
  Butane->AddElement (H, natoms=10);
  
  // Polypropylene
  
  density = 0.9*g/cm3;
  G4Material* Polyprop = new G4Material(name = "Polyprop", density, nel = 2);
  Polyprop->AddElement (C,3);
  Polyprop->AddElement (H,6);

  // Si3N4
  
  density = 3.44*g/cm3;
  G4Material* Si3N4 = new G4Material(name = "Si3N4", density, nel = 2);
  Si3N4->AddElement (Si, natoms=3);
  Si3N4->AddElement (N, natoms=4);
  
  // SiO2
  
  density = 2.5*g/cm3;
  G4Material* SiO2 = new G4Material(name = "SiO2", density, nel = 2);
  SiO2->AddElement (Si, natoms=1);
  SiO2->AddElement (O, natoms=2);
    
  // Brass
  
  density = 8.5*g/cm3;
  G4Material* Laiton = new G4Material(name = "Laiton", density, nel = 2);
  Laiton->AddElement (Cu,1);
  Laiton->AddElement (Zn,1);

  // Phantom
  
  fDensityPhantom = 1.; // in g/cm3

  // Nucleus composition from Alard et al., Rad. Res. 158, 650 (2002) and
  // Comp. Math. Meth. Med. 147252 (2012) 
  //
  // Cytoplasm composition is assumed to be water
  
  // Cytoplasm
  
  fDensityCytoplasm = 1.; // in g/cm3
  density = fDensityCytoplasm*g/cm3;
  G4Material* Cytoplasm1 = new G4Material(name="Cytoplasm1"  , density, ncomponents=2);
  Cytoplasm1->AddElement(H, fractionmass=0.112);
  Cytoplasm1->AddElement(O, fractionmass=0.888);
 
  // Nucleoli
  
  fDensityCytoplasm = 1.;  
  // in g/cm3 (nucleoli are assumed to have the same chemical comp. as nucleus)
  density = fDensityCytoplasm*g/cm3;
  G4Material* Cytoplasm2 = new G4Material(name="Cytoplasm2"  , density, ncomponents=5);
  Cytoplasm2->AddElement(H, fractionmass=0.1064);
  Cytoplasm2->AddElement(O, fractionmass=0.745);
  Cytoplasm2->AddElement(C, fractionmass=0.0904);
  Cytoplasm2->AddElement(N, fractionmass=0.0321);
  Cytoplasm2->AddElement(P, fractionmass=0.0261);
 
  // default is water
  
  fDensityCytoplasm = 1.; // in g/cm3
  density = fDensityCytoplasm*g/cm3;
  G4Material* Cytoplasm3 = new G4Material(name="Cytoplasm3"  , density, ncomponents=2);
  Cytoplasm3->AddElement(H, fractionmass=0.112);
  Cytoplasm3->AddElement(O, fractionmass=0.888);
 
  // Nucleus chemical composition

  fDensityNucleus = 1.; // in g/cm3
  density = fDensityNucleus*g/cm3;
  G4Material* Nucleus1 = new G4Material(name="Nucleus1"  , density, ncomponents=5);
  Nucleus1->AddElement(H, fractionmass=0.1064);
  Nucleus1->AddElement(O, fractionmass=0.745);
  Nucleus1->AddElement(C, fractionmass=0.0904);
  Nucleus1->AddElement(N, fractionmass=0.0321);
  Nucleus1->AddElement(P, fractionmass=0.0261);
 
  fDensityNucleus = 1.; // in g/cm3
  density = fDensityNucleus*g/cm3;
  G4Material* Nucleus2 = new G4Material(name="Nucleus2"  , density, ncomponents=5);
  Nucleus2->AddElement(H, fractionmass=0.1064);
  Nucleus2->AddElement(O, fractionmass=0.745);
  Nucleus2->AddElement(C, fractionmass=0.0904);
  Nucleus2->AddElement(N, fractionmass=0.0321);
  Nucleus2->AddElement(P, fractionmass=0.0261);
 
  // default
  
  fDensityNucleus = 1.; // in g/cm3
  density = fDensityNucleus*g/cm3;
  G4Material* Nucleus3 = new G4Material(name="Nucleus3"  , density, ncomponents=5);
  Nucleus3->AddElement(H, fractionmass=0.1064);
  Nucleus3->AddElement(O, fractionmass=0.745);
  Nucleus3->AddElement(C, fractionmass=0.0904);
  Nucleus3->AddElement(N, fractionmass=0.0321);
  Nucleus3->AddElement(P, fractionmass=0.0261);
 
  // Materials in setup
  
  fDefaultMaterial 	= vacuum;
  fCollimatorMaterial 	= Pt;
  fBoiteMaterial  	= Butane;
  fCathodeMaterial      = Laiton;
  fVerreMaterial        = Si3N4;
  fVerre2Material 	= SiO2;
  fKgmMaterial 		= H2O;
  fBoite2Material       = Air;
  fBoite3Material       = Polyprop;
  
  fNucleusMaterial1 	= Nucleus1;
  fCytoplasmMaterial1 	= Cytoplasm1;
  fNucleusMaterial2 	= Nucleus2;
  fCytoplasmMaterial2 	= Cytoplasm2;
  fNucleusMaterial3 	= Nucleus3;
  fCytoplasmMaterial3 	= Cytoplasm3;
  
  // DISPLAY MATERIALS
  G4cout << G4endl << *(G4Material::GetMaterialTable()) << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* DetectorConstruction::ConstructLine()
{
  // WORLD
  fWorldSizeXY  = 20*m;
  fWorldSizeZ   = 40*m;
   
  // MICROBEAM LINE ANGLE
  fLineAngle = 10*deg;
  
  // TARGET POSITION
  fCiblePositionX = -1461.42*mm;
  fCiblePositionY = 0*mm;
  fCiblePositionZ = -1327 + (955*std::cos(fLineAngle))*mm;
  
  //*************
  // WORLD VOLUME
  //*************
  
  fSolidWorld = new G4Box("World",				         //its name
			  fWorldSizeXY/2,fWorldSizeXY/2,fWorldSizeZ/2);  //its size
  
  
  fLogicWorld = new G4LogicalVolume(fSolidWorld,	//its solid
				    fDefaultMaterial,	//its material
				    "World");		//its name
  
  fPhysiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 "World",		//its name
                                 fLogicWorld,		//its logical volume
                                 NULL,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number

  //*****************
  // FULL LINE VOLUME
  //*****************
  
  fSolidVol = new G4Box("Vol",				      
		        10.*m/2,10.*m/2,(14025)*mm/2);  
		       
  fLogicVol = new G4LogicalVolume(fSolidVol,	       
				  fDefaultMaterial,     
				  "Vol");		
  
  fPhysiVol = new G4PVPlacement(0,			
			       G4ThreeVector(0,0,-2012.5*mm),	
			       "Vol",			
			       fLogicVol,	    
			       fPhysiWorld,	     
			       false,		      
			       0);
			       
  // *************************************************			       		    
  // Whole microbeam line at 10 deg contained in a box 
  // *************************************************

  G4double PosX = fCiblePositionX*mm +( (6958.3/2-3.3)*std::sin(fLineAngle))*mm;
  G4double PosZ = (fCiblePositionZ+2012.5)*mm - ((6958.3/2-3.3)*std::cos(fLineAngle))*mm;

  // Adjust box absolute position
  
  PosX = PosX + 1.3 * micrometer * std::cos(fLineAngle);
  PosZ = PosZ + 1.3 * micrometer * std::sin(fLineAngle);
      
  G4RotationMatrix *rot = new G4RotationMatrix();
  rot->rotateX(0*deg);
  rot->rotateY(10*deg);
  rot->rotateZ(0*deg);
 
  fSolidBoite = new G4Box("Boite", 4*cm, 4*cm, 6958.3*mm/2);
  
  fLogicBoite = new G4LogicalVolume(fSolidBoite, fDefaultMaterial, "Boite");
  
  fPhysiBoite = new G4PVPlacement(rot,
				 G4ThreeVector(PosX,0,PosZ),
				 "Boite", 
				 fLogicBoite,
				 fPhysiVol,
				 false, 
				 0);
  
  //*********************************************************************
  // OBJECT COLLIMATOR (after switching magnet, 5 micrometer in diameter)
  //*********************************************************************
  
  fCollObjSizeXY = 8*cm;
  fCollObjSizeZ = 0.07*mm;
  
  fSolidYoke1 = new G4Box("_CollObj_yoke1_", fCollObjSizeXY/2,fCollObjSizeXY/2,fCollObjSizeZ/2);
  
  fLogicYoke1 = new G4LogicalVolume(fSolidYoke1, fCollimatorMaterial, "_CollObj_yoke1_");
  
  fPhysiYoke1 = new G4PVPlacement( 0, G4ThreeVector(0,0,6958.3*mm/2-3.3*mm-6955*mm+0.07*mm/2), fLogicYoke1, 
                                   "_CollObj_yoke1_",fLogicBoite, false, 0);
   
  // --> FIRST PART
  
  fSolid1Gap = new G4Cons("_CollObj_gap1_", 0.*micrometer, 6*micrometer,
			 0.*micrometer,2.5*micrometer,
			 3.5*micrometer, 
			 0, ((360*CLHEP::pi)/180));
  
  fLogic1Gap = new G4LogicalVolume(fSolid1Gap, fDefaultMaterial, "_CollObj_gap1_");
  
  fPhysi1Gap = new G4PVPlacement(0, G4ThreeVector(0,0,0.0315*mm), fLogic1Gap, "_CollObj_gap1_", 
                                 fLogicYoke1, false, 0);
  
  
  // --> SECOND PART
  
  fSolid2Gap = new G4Cons("_CollObj_gap2_", 0.*micrometer, 15*micrometer,
			 0.*micrometer,6*micrometer,
			 6.5*micrometer, 
			 0, ((360*CLHEP::pi)/180));
  
  fLogic2Gap = new G4LogicalVolume(fSolid2Gap, fDefaultMaterial, "_CollObj_gap2_");
  
  fPhysi2Gap = new G4PVPlacement(0, G4ThreeVector(0,0,0.0215*mm), fLogic2Gap, "_CollObj_gap2_", 
                                 fLogicYoke1, false, 0);
  
  
  // --> THIRD PART
  
  fSolid3Gap = new G4Cons("_CollObj_gap3_", 0.*micrometer, 105*micrometer, 
			 0.*micrometer,15*micrometer,
			 25*micrometer, 
			 0, ((360*CLHEP::pi)/180));
  
  fLogic3Gap = new G4LogicalVolume(fSolid3Gap, fDefaultMaterial, "_CollObj_gap3_");
  
  fPhysi3Gap = new G4PVPlacement(0, G4ThreeVector(0,0,-0.010*mm), fLogic3Gap, "_CollObj_gap3_", fLogicYoke1, 
                                 false, 0);


  //************************
  // GAS DETECTOR COLLIMATOR
  //************************
  
  fSolidYoke2 = new G4Box("_CollDet_yoke_", 2.5*cm, 2.5*cm, 0.035*mm);
 
  fLogicYoke2 = new G4LogicalVolume(fSolidYoke2, fCollimatorMaterial, "_CollDet_yoke_");
  
  fPhysiYoke2 = new G4PVPlacement(0,
                                  G4ThreeVector(0,0,6958.3*mm/2-0.3*mm-3*mm-0.004*mm-0.1*mm-1*mm-2.5*mm-0.070*mm/2), 
                                  fLogicYoke2, "_CollDet_yoke_", fLogicBoite, false, 0);

  // --> FIRST PART
  
  fSolid4Gap = new G4Cons("_CollDet_gap4_", 0.*micrometer, 8*micrometer,
			 0.*micrometer,5*micrometer,
			 7.5*micrometer, 
			 0, ((360*CLHEP::pi)/180));

  fLogic4Gap = new G4LogicalVolume(fSolid4Gap, fDefaultMaterial, "_CollDet_gap4_");
  
  fPhysi4Gap = new G4PVPlacement(0, G4ThreeVector(0,0,0.0275*mm), fLogic4Gap, "_CollDet_gap4_", 
                                 fLogicYoke2, false, 0);
  
  // --> SECOND PART
  
  fSolid5Gap = new G4Cons("_CollDet_gap5_", 0.*micrometer, 105*micrometer,
			 0.*micrometer,8*micrometer,
			 27.5*micrometer, 
			 0, ((360*CLHEP::pi)/180));

  fLogic5Gap = new G4LogicalVolume(fSolid5Gap, fDefaultMaterial, "_CollDet_gap5_");
  
  fPhysi5Gap = new G4PVPlacement(0,
				G4ThreeVector(0,0,-0.0075*mm),
				fLogic5Gap,
				"_CollDet_gap5_", 
				fLogicYoke2,
				false,
				0);
  // ************
  // GAS DETECTOR
  // ************
 
  fSolidBoiteIso = new G4Box("Isobutane", 2.5*cm, 2.5*cm, 1.75*mm);
  
  fLogicBoiteIso = new G4LogicalVolume(fSolidBoiteIso, fBoiteMaterial, "Isobutane");
  
  fPhysiBoiteIso = new G4PVPlacement(0,
				 G4ThreeVector(0,0,6958.3*mm/2-0.3*mm-3*mm-0.004*mm-0.1*mm-3.5*mm/2),
				 "Isobutane", 
				 fLogicBoiteIso,
				 fPhysiBoite,
				 false, 
				 0);
  
  // --> GAS DETECTOR END CAP
  
  fSolidCathode = new G4Box("_Laiton_", 2.5*cm, 2.5*cm, 0.5*mm);
  
  fLogicCathode = new G4LogicalVolume(fSolidCathode, fCathodeMaterial, "_Laiton_");
  
  fPhysiCathode = new G4PVPlacement(0,
				   G4ThreeVector(0,0,1.25*mm),
				   "_Laiton_", 
				   fLogicCathode,
				   fPhysiBoiteIso,
				   false, 0);

  // --> ISOBUTANE GAS  
  
  fSolidIso = new G4Box("_Iso_", 1.*mm, 1.*mm, 0.499925*mm);
  
  fLogicIso = new G4LogicalVolume(fSolidIso, fBoiteMaterial, "_Iso_");
  
  fPhysiIso = new G4PVPlacement(0, 
			       G4ThreeVector(0,0,-0.000075*mm),
			       "_Iso_", 
			       fLogicIso,
			       fPhysiCathode,
			       false, 
			       0);

  // --> Si3N4 WINDOW
  
  fSolidVerre = new G4Box("_Si3N4_", 0.5*mm, 0.5*mm, 0.075*micrometer);
  
  fLogicVerre = new G4LogicalVolume(fSolidVerre, fVerreMaterial, "_Si3N4_");
  
  
  fPhysiVerre = new G4PVPlacement(0,
				 G4ThreeVector(0,0,0.499925*mm),
				 "_Si3N4_", 
				 fLogicVerre,
				 fPhysiCathode,
				 false,
				 0);
  // *******
  // AIR GAP
  // *******
   
  fSolidBoite2 = new G4Box("_Air_", 2.5*cm, 2.5*cm, 0.1*mm/2);
  
  fLogicBoite2 = new G4LogicalVolume(fSolidBoite2, fBoite2Material, "_Air_");
  
  fPhysiBoite2 = new G4PVPlacement(0,
				 G4ThreeVector(0,0,6958.3*mm/2-0.3*mm-3*mm-0.004*mm-0.1*mm/2),
				 "_Air_", 
				 fLogicBoite2,
				 fPhysiBoite,
				 false, 
				 0);

  //*************						
  // CELL SUPPORT
  //*************  
  
  fSolidBoite3 = new G4Box("Polyprop", 2.5*cm, 2.5*cm, 0.004*mm/2);
  
  fLogicBoite3 = new G4LogicalVolume(fSolidBoite3, fBoite3Material, "Polyprop");
  
  fPhysiBoite3 = new G4PVPlacement(0, 
				  G4ThreeVector(0,0,6958.3*mm/2-0.3*mm-3*mm-0.004*mm/2),
				  "Polyprop", 
				  fLogicBoite3, 
				  fPhysiBoite, 
				  false, 
				  0);
  //****
  // KGM   
  //****
    
  fSolidKgm = new G4Box("KGM", 2.5*cm, 2.5*cm, 3*mm/2);

  fLogicKgm = new G4LogicalVolume(fSolidKgm, fKgmMaterial, "KGM");
  
  fPhysiKgm = new G4PVPlacement(0,
			       G4ThreeVector(0,0,6958.3*mm/2-0.3*mm-3*mm/2),
			       "KGM",
			       fLogicKgm,
			       fPhysiBoite, 
			       false,
			       0);

  //*****************
  // MICROSCOPE PLATE
  //*****************
  
  fSolidVerre2 = new G4Box("_Lame_", 2.5*cm, 2.5*cm, 0.150*mm);
  
  fLogicVerre2 = new G4LogicalVolume(fSolidVerre2, fVerre2Material, "_Lame_");
  
  fPhysiVerre2 = new G4PVPlacement(0,
				  G4ThreeVector(0,0,6958.3*mm/2-0.3*mm/2),
				  "_Lame_", 
				  fLogicVerre2,
				  fPhysiBoite, 
				  false,
				  0);

  // **************
  // CELL CYTOPLASM
  // **************
  
  // WITHIN KGM
/*  
  fSolidCyto=new G4Ellipsoid("CYTO",25*micrometer, 25*micrometer, 11*micrometer);
 
  fLogicCyto=new G4LogicalVolume (fSolidCyto, fDefaultMaterial, "CYTO");

  fPhysiCyto=new G4PVPlacement(0, G4ThreeVector(0,0,-1.5*mm+11*micrometer),"CYTO",fLogicCyto, fPhysiKgm, false, 0);
*/

  // ************
  // CELL PHANTOM
  // ************

  fMyCellParameterisation = new CellParameterisation
        (fNucleusMaterial1,fCytoplasmMaterial1,
	 fNucleusMaterial2,fCytoplasmMaterial2,
	 fNucleusMaterial3,fCytoplasmMaterial3);

  fSolidPhantom = new G4Box("Phantom", 
  	fMyCellParameterisation->GetPixelSizeX()/2, 
	fMyCellParameterisation->GetPixelSizeY()/2, 
	fMyCellParameterisation->GetPixelSizeZ()/2); 
  
  fLogicPhantom = new G4LogicalVolume(fSolidPhantom,fDefaultMaterial,"Phantom",0,0,0);
    
  SetNbOfPixelsInPhantom (fMyCellParameterisation->GetPhantomTotalPixels());

  SetMassNucleus(fMyCellParameterisation->GetNucleusMass());

  SetMassCytoplasm(fMyCellParameterisation->GetCytoplasmMass());

  fPhysiPhantom = new G4PVParameterised(
                            "Phantom",        // their name
                            fLogicPhantom,    // their logical volumr
                            //logicCyto,      // Mother logical volume is Cyto
                            fLogicKgm,        // Mother logical volume is Kgm
			    kUndefined,       // Are placed along this axis 
                            fMyCellParameterisation->GetPhantomTotalPixels(),    // Number of boxes
                            fMyCellParameterisation,false);   // The parametrisation

  G4cout << " ==========> The phantom contains " << fMyCellParameterisation->GetPhantomTotalPixels() << " voxels " << G4endl;		    
  G4cout << " ==========> Nucleus mass (kg)=" << fMyCellParameterisation->GetNucleusMass() / kg << G4endl;
  G4cout << " ==========> Cytoplasm mass (kg)=" << fMyCellParameterisation->GetCytoplasmMass()/ kg << G4endl;
  G4cout << " ==========> Voxel size X (um)=" << fMyCellParameterisation->GetPixelSizeX()/um << G4endl;
  G4cout << " ==========> Voxel size Y (um)=" << fMyCellParameterisation->GetPixelSizeY()/um << G4endl;
  G4cout << " ==========> Voxel size Z (um)=" << fMyCellParameterisation->GetPixelSizeZ()/um << G4endl; 
  G4cout << G4endl; 
				    		    
  // USER LIMITS ON STEP LENGTH
  
/*
  fLogicWorld->SetUserLimits(new G4UserLimits(100*mm));
  fLogicVol->SetUserLimits(new G4UserLimits(100*mm));
  fLogicBoite->SetUserLimits(new G4UserLimits(10*mm));
*/

/*
  logicPhantom->SetUserLimits (new G4UserLimits(0.5*micrometer));
  logic1Gap->SetUserLimits (new G4UserLimits(5*micrometer));
  logic2Gap->SetUserLimits (new G4UserLimits(5*micrometer));
  logic3Gap->SetUserLimits (new G4UserLimits(5*micrometer));
  logic4Gap->SetUserLimits (new G4UserLimits(5*micrometer));
  logic5Gap->SetUserLimits (new G4UserLimits(5*micrometer));
  logicBoiteIso->SetUserLimits (new G4UserLimits(200.*micrometer));
  logicCathode->SetUserLimits (new G4UserLimits(100.*micrometer));
  logicIso->SetUserLimits (new G4UserLimits(100.*micrometer));
  logicVerre->SetUserLimits (new G4UserLimits(0.02*micrometer));
  logicBoite2->SetUserLimits (new G4UserLimits(10*micrometer));
  logicBoite3->SetUserLimits (new G4UserLimits(0.2*micrometer));
  logicKgm->SetUserLimits (new G4UserLimits(1*micrometer));
  logicVerre2->SetUserLimits (new G4UserLimits(10*micrometer));
*/

  // Relaxed 
  fLogicWorld->SetUserLimits(new G4UserLimits(10*mm));
  fLogicVol->SetUserLimits(new G4UserLimits(10*mm));
  fLogicBoite->SetUserLimits(new G4UserLimits(1*mm));

  // VISUALISATION ATTRIBUTES (for phantom, see in Parameterisation class)
  
  G4VisAttributes* simpleWorldVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0)); //White
  simpleWorldVisAtt->SetVisibility(true);
  
  G4VisAttributes* simplePlain= new G4VisAttributes(G4Colour(1.0,1.0,1.0)); //White
  simplePlain->SetVisibility(true);
  simplePlain->SetForceSolid(true);
    
  G4VisAttributes* simpleBoxAttLine= new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  simpleBoxAttLine->SetVisibility(true);
   
  G4VisAttributes* simpleBoxAtt= new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  simpleBoxAtt->SetDaughtersInvisible(false);
  simpleBoxAtt->SetForceSolid(false);
 
  G4VisAttributes* simpleBoxAtt2= new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  simpleBoxAtt2->SetDaughtersInvisible(false);
  simpleBoxAtt2->SetForceSolid(false);
  
  G4VisAttributes* simpleBoxAttKGM= new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  simpleBoxAttKGM->SetDaughtersInvisible(false);
  simpleBoxAttKGM->SetForceSolid(false);
  
  G4VisAttributes* simpleBoxAttPropyl= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxAttPropyl->SetDaughtersInvisible(true);
  simpleBoxAttPropyl->SetForceSolid(false);
  
  G4VisAttributes* simpleBoxAttAir= new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  simpleBoxAttAir->SetDaughtersInvisible(true);
  simpleBoxAttAir->SetForceSolid(false);
  
  G4VisAttributes* simpleBoxAtt3= new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  simpleBoxAtt3->SetDaughtersInvisible(false);
  simpleBoxAtt3->SetForceSolid(false);
  
  fLogicYoke1->SetVisAttributes(simpleBoxAtt);
  fLogic1Gap->SetVisAttributes(simpleBoxAtt);
  fLogic2Gap->SetVisAttributes(simpleBoxAtt);
  fLogic3Gap->SetVisAttributes(simpleBoxAtt);
  fLogicYoke2->SetVisAttributes(simpleBoxAtt);  
  fLogic4Gap->SetVisAttributes(simpleBoxAtt);
  fLogic5Gap->SetVisAttributes(simpleBoxAtt);
  fLogicBoite->SetVisAttributes(simpleBoxAttLine);
  fLogicCathode->SetVisAttributes(simpleBoxAttPropyl);
  fLogicIso->SetVisAttributes(simpleBoxAttPropyl);
  fLogicBoiteIso->SetVisAttributes(simpleBoxAttPropyl);
  fLogicVerre->SetVisAttributes(simpleBoxAtt);
  fLogicBoite2->SetVisAttributes(simpleBoxAttAir);
  fLogicBoite3->SetVisAttributes(simpleBoxAtt);
  fLogicKgm->SetVisAttributes(simpleBoxAttKGM);
  fLogicVerre2->SetVisAttributes(simpleBoxAtt);
  
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::ConstructSDandField()
{
  if(!fField) fField = new EMField(); 
  
  G4EqMagElectricField* fEquation = new G4EqMagElectricField(fField);
  G4MagIntegratorStepper* fStepper = new G4ClassicalRK4 (fEquation,8);
  G4FieldManager* fFieldMgr = 
    G4TransportationManager::GetTransportationManager()->GetFieldManager();

  // Relaxed
  G4MagInt_Driver* fIntgrDriver = 
    new G4MagInt_Driver(1*mm,fStepper,fStepper->GetNumberOfVariables() );

  G4ChordFinder* fChordFinder = new G4ChordFinder(fIntgrDriver);
  fFieldMgr->SetChordFinder(fChordFinder);
  fFieldMgr->SetDetectorField(fField);

  // FOLLOWING PARAMETERS TUNED FROM RAY-TRACING SIMULATIONS OF THE AIFIRA NANOBEAM LINE
  /*
  fFieldMgr->GetChordFinder()->SetDeltaChord(1e-9*m);
  fFieldMgr->SetDeltaIntersection(1e-9*m);
  fFieldMgr->SetDeltaOneStep(1e-9*m);     
      
  fPropInField =
    G4TransportationManager::GetTransportationManager()->GetPropagatorInField();
  fPropInField->SetMinimumEpsilonStep(1e-16); // instead of 11
  fPropInField->SetMaximumEpsilonStep(1e-15); // instead of 10
  */
}
