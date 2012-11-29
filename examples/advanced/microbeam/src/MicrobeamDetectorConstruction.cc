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
// -------------------------------------------------------------------
// $Id$
// -------------------------------------------------------------------

#include "MicrobeamDetectorConstruction.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrobeamDetectorConstruction::MicrobeamDetectorConstruction()
  
  :defaultMaterial(NULL),collimatorMaterial(NULL),BoiteMaterial(NULL),
   CathodeMaterial(NULL),VerreMaterial(NULL),Verre2Material(NULL),   
   KgmMaterial(NULL),Boite2Material(NULL),Boite3Material(NULL),
   nucleusMaterial1(NULL),cytoplasmMaterial1(NULL),
   nucleusMaterial2(NULL),cytoplasmMaterial2(NULL),
   nucleusMaterial3(NULL),cytoplasmMaterial3(NULL),
   physiWorld(NULL),logicWorld(NULL),solidWorld(NULL),
   physiVol(NULL),logicVol(NULL),solidVol(NULL),
   physiBoite(NULL),logicBoite(NULL),solidBoite(NULL),
   physiYoke1(NULL),logicYoke1(NULL),solidYoke1(NULL),
   physi1Gap(NULL),logic1Gap(NULL),solid1Gap(NULL),
   physi2Gap(NULL),logic2Gap(NULL),solid2Gap(NULL), 
   physi3Gap(NULL),logic3Gap(NULL),solid3Gap(NULL),
   physiYoke2(NULL),logicYoke2(NULL),solidYoke2(NULL),
   physi4Gap(NULL),logic4Gap(NULL),solid4Gap(NULL),
   physi5Gap(NULL),logic5Gap(NULL),solid5Gap(NULL), 
   physiBoiteIso(NULL),logicBoiteIso(NULL),solidBoiteIso(NULL),
   physiCathode(NULL),logicCathode(NULL),solidCathode(NULL), 
   physiIso(NULL),logicIso(NULL),solidIso(NULL),
   physiVerre(NULL),logicVerre(NULL),solidVerre(NULL),
   physiBoite2(NULL),logicBoite2(NULL),solidBoite2(NULL),
   physiBoite3(NULL),logicBoite3(NULL),solidBoite3(NULL),
   physiKgm(NULL),logicKgm(NULL),solidKgm(NULL),
   physiVerre2(NULL),logicVerre2(NULL),solidVerre2(NULL),
   physiPhantom(NULL),logicPhantom(NULL),solidPhantom(NULL)
  
{
  WorldSizeXY=WorldSizeZ=0;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrobeamDetectorConstruction::~MicrobeamDetectorConstruction()
{
  delete defaultMaterial;
  delete collimatorMaterial;
  delete BoiteMaterial;
  delete CathodeMaterial;
  delete VerreMaterial;
  delete Verre2Material;
  delete KgmMaterial;
  delete Boite2Material;
  delete Boite3Material;
  delete nucleusMaterial1;
  delete cytoplasmMaterial1;
  delete nucleusMaterial2;
  delete cytoplasmMaterial2;
  delete nucleusMaterial3;
  delete cytoplasmMaterial3;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* MicrobeamDetectorConstruction::Construct()
  
{
  DefineMaterials();
  return ConstructMicrobeamLine();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicrobeamDetectorConstruction::DefineMaterials()
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
 
  // Vaccum standard definition...
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
  
  // Low Pressure Air
  density = (5e-6/1013.)*1.290*mg/cm3; // 5e-6 mbar is the usual beam pipe air pressure
  pressure = 1*atmosphere;
  temperature = 293.16*kelvin;
  G4Material* LPAir = new G4Material(name="LPAir"  , density, ncomponents=3, kStateGas, temperature, pressure);
  LPAir->AddElement(N, fractionmass=0.715);
  LPAir->AddElement(O, fractionmass=0.25);
  LPAir->AddElement(Ar, fractionmass=0.035);
  
  // Platine
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
    
  // Laiton
  density = 8.5*g/cm3;
  G4Material* Laiton = new G4Material(name = "Laiton", density, nel = 2);
  Laiton->AddElement (Cu,1);
  Laiton->AddElement (Zn,1);

  // Phantom
  densityPhantom = 1.; // in g/cm3

  // Nucleus composition from Alard et al., Rad. Res. 158, 650 (2002)

  // Cytoplasm chemical composition
  densityCytoplasm = 1.; // in g/cm3
  density = densityCytoplasm*g/cm3;
  G4Material* Cytoplasm1 = new G4Material(name="Cytoplasm1"  , density, ncomponents=2);
  Cytoplasm1->AddElement(H, fractionmass=0.112);
  Cytoplasm1->AddElement(O, fractionmass=0.888);
 
  densityCytoplasm = 1.; 
  // in g/cm3 (nucleoli are assumed to have the same chemical comp. as nucleus)
  density = densityCytoplasm*g/cm3;
  G4Material* Cytoplasm2 = new G4Material(name="Cytoplasm2"  , density, ncomponents=5);
  Cytoplasm2->AddElement(H, fractionmass=0.1064);
  Cytoplasm2->AddElement(O, fractionmass=0.745);
  Cytoplasm2->AddElement(C, fractionmass=0.0904);
  Cytoplasm2->AddElement(N, fractionmass=0.0321);
  Cytoplasm2->AddElement(P, fractionmass=0.0261);
 
  // default 
  densityCytoplasm = 1.; // in g/cm3
  density = densityCytoplasm*g/cm3;
  G4Material* Cytoplasm3 = new G4Material(name="Cytoplasm3"  , density, ncomponents=2);
  Cytoplasm3->AddElement(H, fractionmass=0.112);
  Cytoplasm3->AddElement(O, fractionmass=0.888);
 
  // Nucleus chemical composition
  densityNucleus = 1.; // in g/cm3
  density = densityNucleus*g/cm3;
  G4Material* Nucleus1 = new G4Material(name="Nucleus1"  , density, ncomponents=5);
  Nucleus1->AddElement(H, fractionmass=0.1064);
  Nucleus1->AddElement(O, fractionmass=0.745);
  Nucleus1->AddElement(C, fractionmass=0.0904);
  Nucleus1->AddElement(N, fractionmass=0.0321);
  Nucleus1->AddElement(P, fractionmass=0.0261);
 
  densityNucleus = 1.; // in g/cm3
  density = densityNucleus*g/cm3;
  G4Material* Nucleus2 = new G4Material(name="Nucleus2"  , density, ncomponents=5);
  Nucleus2->AddElement(H, fractionmass=0.1064);
  Nucleus2->AddElement(O, fractionmass=0.745);
  Nucleus2->AddElement(C, fractionmass=0.0904);
  Nucleus2->AddElement(N, fractionmass=0.0321);
  Nucleus2->AddElement(P, fractionmass=0.0261);
 
  // default
  densityNucleus = 1.; // in g/cm3
  density = densityNucleus*g/cm3;
  G4Material* Nucleus3 = new G4Material(name="Nucleus3"  , density, ncomponents=5);
  Nucleus3->AddElement(H, fractionmass=0.1064);
  Nucleus3->AddElement(O, fractionmass=0.745);
  Nucleus3->AddElement(C, fractionmass=0.0904);
  Nucleus3->AddElement(N, fractionmass=0.0321);
  Nucleus3->AddElement(P, fractionmass=0.0261);
 
  // Materials in setup.
  defaultMaterial 	= vacuum;
  collimatorMaterial 	= Pt;
  BoiteMaterial  	= Butane;
  CathodeMaterial       = Laiton;
  VerreMaterial         = Si3N4;
  Verre2Material 	= SiO2;
  KgmMaterial 		= H2O;
  Boite2Material        = Air;
  Boite3Material        = Polyprop;
  nucleusMaterial1 	= Nucleus1;
  cytoplasmMaterial1 	= Cytoplasm1;
  nucleusMaterial2 	= Nucleus2;
  cytoplasmMaterial2 	= Cytoplasm2;
  nucleusMaterial3 	= Nucleus3;
  cytoplasmMaterial3 	= Cytoplasm3;
  
  // DISPLAY MATERIALS
  G4cout << G4endl << *(G4Material::GetMaterialTable()) << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* MicrobeamDetectorConstruction::ConstructMicrobeamLine()
{
  // WORLD
  WorldSizeXY  = 20*m;
  WorldSizeZ   = 40*m;
   
  // MICROBEAM LINE ANGLE
  lineAngle = 10*deg;
  
  // TARGET POSITION
  CiblePositionX = -1461.42*mm;
  CiblePositionY = 0*mm;
  CiblePositionZ = -1327 + (955*std::cos(lineAngle))*mm;
  
  // ELECTROMAGNETIC FIELD PARAMETERS

  static G4bool fieldIsInitialized = false;
  
  if(!fieldIsInitialized)

    {
      Field = new MicrobeamEMField();
      pEquation = new G4EqMagElectricField(Field);
      pStepper = new G4ClassicalRK4 (pEquation);
      pFieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
      pIntgrDriver = new G4MagInt_Driver(0.000001*mm,pStepper,pStepper->GetNumberOfVariables() );
      pChordFinder = new G4ChordFinder(pIntgrDriver);
      pFieldMgr->SetChordFinder( pChordFinder );
      pFieldMgr->SetDetectorField(Field);
      fieldIsInitialized = true;
      
      // FOLLOWING PARAMETERS TUNED FROM RAY-TRACING SIMULATIONS OF THE AIFIRA NANOBEAM LINE
      
      pFieldMgr->GetChordFinder()->SetDeltaChord(1e-9*m);
      pFieldMgr->SetDeltaIntersection(1e-9*m);
      pFieldMgr->SetDeltaOneStep(1e-9*m);     
      
      propInField =
	G4TransportationManager::GetTransportationManager()->GetPropagatorInField();
      propInField->SetMinimumEpsilonStep(1e-11);
      propInField->SetMaximumEpsilonStep(1e-10);
    }

  //*************
  // WORLD VOLUME
  //*************
  
  solidWorld = new G4Box("World",				       //its name
			 WorldSizeXY/2,WorldSizeXY/2,WorldSizeZ/2);  //its size
  
  
  logicWorld = new G4LogicalVolume(solidWorld,	        //its solid
				   defaultMaterial,	//its material
				   "World");		//its name
  
  physiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 "World",		//its name
                                 logicWorld,		//its logical volume
                                 NULL,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number

  //*****************
  // FULL LINE VOLUME
  //*****************
  
  solidVol = new G4Box("Vol",				      
		       10.*m/2,10.*m/2,(14025)*mm/2);  
		       
  logicVol = new G4LogicalVolume(solidVol,	       
				 defaultMaterial,     
				 "Vol");		
  
  physiVol = new G4PVPlacement(0,			
			       G4ThreeVector(0,0,-2012.5*mm),	
			       "Vol",			
			       logicVol,	    
			       physiWorld,	     
			       false,		      
			       0);
			       
  // *************************************************			       		    
  // Whole microbeam line at 10 deg contained in a box 
  // *************************************************

  G4double PosX = CiblePositionX*mm +( (6958.3/2-3.3)*std::sin(lineAngle))*mm;
  G4double PosZ = (CiblePositionZ+2012.5)*mm - ((6958.3/2-3.3)*std::cos(lineAngle))*mm;

  // Adjust box absolute position
  
  PosX = PosX + 1.3 * micrometer * std::cos(lineAngle);
  PosZ = PosZ + 1.3 * micrometer * std::sin(lineAngle);
      
  G4RotationMatrix *rot = new G4RotationMatrix();
  rot->rotateX(0*deg);
  rot->rotateY(10*deg);
  rot->rotateZ(0*deg);
 
  solidBoite = new G4Box("Boite", 4*cm, 4*cm, 6958.3*mm/2);
  
  logicBoite = new G4LogicalVolume(solidBoite, defaultMaterial, "Boite");
  
  physiBoite = new G4PVPlacement(rot,
				 G4ThreeVector(PosX,0,PosZ),
				 "Boite", 
				 logicBoite,
				 physiVol,
				 false, 
				 0);
  
  //*********************************************************************
  // OBJECT COLLIMATOR (after switching magnet, 5 micrometer in diameter)
  //*********************************************************************
  
  CollObjSizeXY = 8*cm;
  CollObjSizeZ = 0.07*mm;
  
  solidYoke1 = new G4Box("_CollObj_yoke1_", CollObjSizeXY/2,CollObjSizeXY/2 , CollObjSizeZ/2);
  
  logicYoke1 = new G4LogicalVolume(solidYoke1, collimatorMaterial, "_CollObj_yoke1_");
  physiYoke1 = new G4PVPlacement( 0, G4ThreeVector(0,0,6958.3*mm/2-3.3*mm-6955*mm+0.07*mm/2), logicYoke1, "_CollObj_yoke1_", 
				  logicBoite, false, 0);
   
  // --> FIRST PART
  
  solid1Gap = new G4Cons("_CollObj_gap1_", 0.*micrometer, 6*micrometer,
			 0.*micrometer,2.5*micrometer,
			 3.5*micrometer, 
			 0, ((360*CLHEP::pi)/180));
  
  logic1Gap = new G4LogicalVolume(solid1Gap, defaultMaterial, "_CollObj_gap1_");
  
  physi1Gap = new G4PVPlacement(0, G4ThreeVector(0,0,0.0315*mm), logic1Gap, "_CollObj_gap1_", logicYoke1, false, 0);
  
  
  // --> SECOND PART
  
  solid2Gap = new G4Cons("_CollObj_gap2_", 0.*micrometer, 15*micrometer,
			 0.*micrometer,6*micrometer,
			 6.5*micrometer, 
			 0, ((360*CLHEP::pi)/180));
  
  logic2Gap = new G4LogicalVolume(solid2Gap, defaultMaterial, "_CollObj_gap2_");
  
  physi2Gap = new G4PVPlacement(0, G4ThreeVector(0,0,0.0215*mm), logic2Gap, "_CollObj_gap2_", logicYoke1, false, 0);
  
  
  // --> THIRD PART
  
  solid3Gap = new G4Cons("_CollObj_gap3_", 0.*micrometer, 105*micrometer, 
			 0.*micrometer,15*micrometer,
			 25*micrometer, 
			 0, ((360*CLHEP::pi)/180));
  
  logic3Gap = new G4LogicalVolume(solid3Gap, defaultMaterial, "_CollObj_gap3_");
  
  physi3Gap = new G4PVPlacement(0, G4ThreeVector(0,0,-0.010*mm), logic3Gap, "_CollObj_gap3_", logicYoke1, false, 0);


  //************************
  // GAS DETECTOR COLLIMATOR
  //************************
  
  solidYoke2 = new G4Box("_CollDet_yoke_", 2.5*cm, 2.5*cm, 0.035*mm);
 
  logicYoke2 = new G4LogicalVolume(solidYoke2, collimatorMaterial, "_CollDet_yoke_");
  
  physiYoke2 = new G4PVPlacement(0, G4ThreeVector(0,0,6958.3*mm/2-0.3*mm-3*mm-0.004*mm-0.1*mm-1*mm-2.5*mm-0.070*mm/2), logicYoke2, "_CollDet_yoke_", 
				 logicBoite, false, 0);

  // --> FIRST PART
  
  solid4Gap = new G4Cons("_CollDet_gap4_", 0.*micrometer, 8*micrometer,
			 0.*micrometer,5*micrometer,
			 7.5*micrometer, 
			 0, ((360*CLHEP::pi)/180));

  logic4Gap = new G4LogicalVolume(solid4Gap, defaultMaterial, "_CollDet_gap4_");
  
  physi4Gap = new G4PVPlacement(0, G4ThreeVector(0,0,0.0275*mm), logic4Gap, "_CollDet_gap4_", logicYoke2, false, 0);
  
  // --> SECOND PART
  
  solid5Gap = new G4Cons("_CollDet_gap5_", 0.*micrometer, 105*micrometer,
			 0.*micrometer,8*micrometer,
			 27.5*micrometer, 
			 0, ((360*CLHEP::pi)/180));

  logic5Gap = new G4LogicalVolume(solid5Gap, defaultMaterial, "_CollDet_gap5_");
  
  physi5Gap = new G4PVPlacement(0,
				G4ThreeVector(0,0,-0.0075*mm),
				logic5Gap,
				"_CollDet_gap5_", 
				logicYoke2,
				false,
				0);
  // ************
  // GAS DETECTOR
  // ************
 
  solidBoiteIso = new G4Box("Isobutane", 2.5*cm, 2.5*cm, 1.75*mm);
  
  logicBoiteIso = new G4LogicalVolume(solidBoiteIso, BoiteMaterial, "Isobutane");
  
  physiBoiteIso = new G4PVPlacement(0,
				 G4ThreeVector(0,0,6958.3*mm/2-0.3*mm-3*mm-0.004*mm-0.1*mm-3.5*mm/2),
				 "Isobutane", 
				 logicBoiteIso,
				 physiBoite,
				 false, 
				 0);
  // --> GAS DETECTOR END CAP
  
  solidCathode = new G4Box("_Laiton_", 2.5*cm, 2.5*cm, 0.5*mm);
  
  logicCathode = new G4LogicalVolume(solidCathode, CathodeMaterial, "_Laiton_");
  
  physiCathode = new G4PVPlacement(0,
				   G4ThreeVector(0,0,1.25*mm),
				   "_Laiton_", 
				   logicCathode,
				   physiBoiteIso,
				   false, 0);

  // --> ISOBUTANE GAS  
  
  solidIso = new G4Box("_Iso_", 1.*mm, 1.*mm, 0.499925*mm);
  
  logicIso = new G4LogicalVolume(solidIso, BoiteMaterial, "_Iso_");
  
  physiIso = new G4PVPlacement(0, 
			       G4ThreeVector(0,0,-0.000075*mm),
			       "_Iso_", 
			       logicIso,
			       physiCathode,
			       false, 
			       0);

  // --> Si3N4 WINDOW
  
  solidVerre = new G4Box("_Si3N4_", 0.5*mm, 0.5*mm, 0.075*micrometer);
  
  logicVerre = new G4LogicalVolume(solidVerre, VerreMaterial, "_Si3N4_");
  
  
   physiVerre = new G4PVPlacement(0,
				 G4ThreeVector(0,0,0.499925*mm),
				 "_Si3N4_", 
				 logicVerre,
				 physiCathode,
				 false,
				 0);
  // *******
  // AIR GAP
  // *******
   
  solidBoite2 = new G4Box("_Air_", 2.5*cm, 2.5*cm, 0.1*mm/2);
  
  logicBoite2 = new G4LogicalVolume(solidBoite2, Boite2Material, "_Air_");
  
  physiBoite2 = new G4PVPlacement(0,
				 G4ThreeVector(0,0,6958.3*mm/2-0.3*mm-3*mm-0.004*mm-0.1*mm/2),
				 "_Air_", 
				 logicBoite2,
				 physiBoite,
				 false, 
				 0);

  //*************						
  // CELL SUPPORT
  //*************  
  
  solidBoite3 = new G4Box("Polyprop", 2.5*cm, 2.5*cm, 0.004*mm/2);
  
  logicBoite3 = new G4LogicalVolume(solidBoite3, Boite3Material, "Polyprop");
  
  physiBoite3 = new G4PVPlacement(0, 
				  G4ThreeVector(0,0,6958.3*mm/2-0.3*mm-3*mm-0.004*mm/2),
				  "Polyprop", 
				  logicBoite3, 
				  physiBoite, 
				  false, 
				  0);
  //****
  // KGM   
  //****
    
  solidKgm = new G4Box("KGM", 2.5*cm, 2.5*cm, 3*mm/2);

  logicKgm = new G4LogicalVolume(solidKgm, KgmMaterial, "KGM");
  
  physiKgm = new G4PVPlacement(0,
			       G4ThreeVector(0,0,6958.3*mm/2-0.3*mm-3*mm/2),
			       "KGM",
			       logicKgm,
			       physiBoite, 
			       false,
			       0);

  //*****************
  // MICROSCOPE PLATE
  //*****************
  
  solidVerre2 = new G4Box("_Lame_", 2.5*cm, 2.5*cm, 0.150*mm);
  
  logicVerre2 = new G4LogicalVolume(solidVerre2, Verre2Material, "_Lame_");
  
  physiVerre2 = new G4PVPlacement(0,
				  G4ThreeVector(0,0,6958.3*mm/2-0.3*mm/2),
				  "_Lame_", 
				  logicVerre2,
				  physiBoite, 
				  false,
				  0);

  // **************
  // CELL CYTOPLASM
  // **************
  
  // WITHIN KGM
/*  
  solidCyto=new G4Ellipsoid("CYTO",25*micrometer, 25*micrometer, 11*micrometer);
 
  logicCyto=new G4LogicalVolume (solidCyto, defaultMaterial, "CYTO");

  physiCyto=new G4PVPlacement(0, G4ThreeVector(0,0,-1.5*mm+11*micrometer),"CYTO",logicCyto,physiKgm, false, 0);
*/

  // ************
  // CELL PHANTOM
  // ************

  solidPhantom = new G4Box("Phantom", 
  	myMicrobeamPhantomConfiguration.GetPixelSizeX()/2, 
	myMicrobeamPhantomConfiguration.GetPixelSizeY()/2, 
	myMicrobeamPhantomConfiguration.GetPixelSizeZ()/2); 
  
  logicPhantom = new G4LogicalVolume(solidPhantom,defaultMaterial,"Phantom",0,0,0);
    
    // PHANTOM MASSES

  SetNbOfPixelsInPhantom (myMicrobeamPhantomConfiguration.GetPhantomTotalPixels());

  SetMassNucleus(myMicrobeamPhantomConfiguration.GetNucleusMass());

  SetMassCytoplasm(myMicrobeamPhantomConfiguration.GetCytoplasmMass());

  // PHANTOM

  phantomParam = new MicrobeamCellParameterisation
  	(myMicrobeamPhantomConfiguration.GetPhantomTotalPixels(),
	 myMicrobeamPhantomConfiguration.GetPixelSizeX()/2,
	 myMicrobeamPhantomConfiguration.GetPixelSizeY()/2,
	 myMicrobeamPhantomConfiguration.GetPixelSizeZ()/2,
	 nucleusMaterial1,cytoplasmMaterial1,
	 nucleusMaterial2,cytoplasmMaterial2,
	 nucleusMaterial3,cytoplasmMaterial3
	 );

  physiPhantom = new G4PVParameterised(
                            "Phantom",       // their name
                            logicPhantom,    // their logical volumr
//                            logicCyto,       // Mother logical volume
                            logicKgm,       // Mother logical volume
			    kUndefined,          // Are placed along this axis 
                            phantomParam->GetNoBoxes(),    // Number of boxes
                            phantomParam,false);   // The parametrisation

  G4cout << " ==========> The phantom contains " 
    << myMicrobeamPhantomConfiguration.GetPhantomTotalPixels() << " voxels " << G4endl;		    
  G4cout << G4endl; 
				    		    
  // USER LIMITS ON STEP LENGTH
  
  logicWorld->SetUserLimits(new G4UserLimits(100*mm));
  logicVol->SetUserLimits(new G4UserLimits(100*mm));
  logicBoite->SetUserLimits(new G4UserLimits(10*mm));

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
  
  logicYoke1->SetVisAttributes(simpleBoxAtt);
  logic1Gap->SetVisAttributes(simpleBoxAtt);
  logic2Gap->SetVisAttributes(simpleBoxAtt);
  logic3Gap->SetVisAttributes(simpleBoxAtt);
  logicYoke2->SetVisAttributes(simpleBoxAtt);  
  logic4Gap->SetVisAttributes(simpleBoxAtt);
  logic5Gap->SetVisAttributes(simpleBoxAtt);
  logicBoite->SetVisAttributes(simpleBoxAttLine);
  logicCathode->SetVisAttributes(simpleBoxAttPropyl);
  logicIso->SetVisAttributes(simpleBoxAttPropyl);
  logicBoiteIso->SetVisAttributes(simpleBoxAttPropyl);
  logicVerre->SetVisAttributes(simpleBoxAtt);
  logicBoite2->SetVisAttributes(simpleBoxAttAir);
  logicBoite3->SetVisAttributes(simpleBoxAtt);
  logicKgm->SetVisAttributes(simpleBoxAttKGM);
  logicVerre2->SetVisAttributes(simpleBoxAtt);

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
  return physiWorld;
}
