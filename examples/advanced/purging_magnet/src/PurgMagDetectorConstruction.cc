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
// Code developed by:
//  S.Larsson
//
//    *****************************************
//    *                                       *
//    *    PurgMagDetectorConstruction.cc     *
//    *                                       *
//    *****************************************
//
//
#include "PurgMagDetectorConstruction.hh"
#include "PurgMagTabulatedField3D.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVParameterised.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EqMagElectricField.hh"

#include "G4ChordFinder.hh"
#include "G4UniformMagField.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Possibility to turn off (0) magnetic field and measurement volume. 
#define GAP 1          // Magnet geometric volume
#define MAG 1          // Magnetic field grid
#define MEASUREVOL 1   // Volume for measurement

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PurgMagDetectorConstruction::PurgMagDetectorConstruction()

  :physiWorld(NULL), logicWorld(NULL), solidWorld(NULL),
   physiGap1(NULL), logicGap1(NULL), solidGap1(NULL),
   physiGap2(NULL), logicGap2(NULL), solidGap2(NULL),
   physiMeasureVolume(NULL), logicMeasureVolume(NULL), 
   solidMeasureVolume(NULL),
   WorldMaterial(NULL), 
   GapMaterial(NULL)
    
{
  fField.Put(0);
  WorldSizeXY=WorldSizeZ=0;
  GapSizeX1=GapSizeX2=GapSizeY1=GapSizeY2=GapSizeZ=0;
  MeasureVolumeSizeXY=MeasureVolumeSizeZ=0;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PurgMagDetectorConstruction::~PurgMagDetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* PurgMagDetectorConstruction::Construct()

{
  DefineMaterials();
  return ConstructCalorimeter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PurgMagDetectorConstruction::DefineMaterials()
{ 
  //This function illustrates the possible ways to define materials.
  //Density and mass per mole taken from Physics Handbook for Science
  //and engineering, sixth edition. This is a general material list
  //with extra materials for other examples.
  
  G4String name, symbol;             
  G4double density;            
  
  G4int ncomponents, natoms;
  G4double fractionmass;
  G4double temperature, pressure;
  
  // Define Elements  
  // Example: G4Element* Notation  = new G4Element ("Element", "Notation", z, a);
  G4Element*   H  = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);
  G4Element*   N  = new G4Element ("Nitrogen", "N", 7., 14.01*g/mole);
  G4Element*   O  = new G4Element ("Oxygen"  , "O", 8. , 16.00*g/mole);
  G4Element*   Ar = new G4Element ("Argon" , "Ar", 18., 39.948*g/mole );
  
  
  // Define Material
  // Example: G4Material* Notation = new G4Material("Material", z, a, density);
  /* Not used in this setup, will be used in further development.
  G4Material* He = new G4Material("Helium", 2., 4.00*g/mole, 0.178*mg/cm3);
  G4Material* Be = new G4Material("Beryllium", 4., 9.01*g/mole, 1.848*g/cm3);
  G4Material* W  = new G4Material("Tungsten", 74., 183.85*g/mole, 19.30*g/cm3);
  G4Material* Cu = new G4Material("Copper", 29., 63.55*g/mole, 8.96*g/cm3);
  */
  G4Material* Fe = new G4Material("Iron", 26., 55.84*g/mole, 7.87*g/cm3);  

  // Define materials from elements.
  
  // Case 1: chemical molecule  
  // Water 
  density = 1.000*g/cm3;
  G4Material* H2O = new G4Material(name="H2O"  , density, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  
  // Case 2: mixture by fractional mass.
  // Air
  density = 1.290*mg/cm3;
  G4Material* Air = new G4Material(name="Air"  , density, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  // Vacuum
  density     = 1.e-5*g/cm3;
  pressure    = 2.e-2*bar;
  temperature = STP_Temperature;         //from PhysicalConstants.h
  G4Material* vacuum = new G4Material(name="vacuum", density, ncomponents=1,
                                      kStateGas,temperature,pressure);
  vacuum->AddMaterial(Air, fractionmass=1.);


  // Laboratory vacuum: Dry air (average composition)
  density = 1.7836*mg/cm3 ;       // STP
  G4Material* Argon = new G4Material(name="Argon", density, ncomponents=1);
  Argon->AddElement(Ar, 1);
  
  density = 1.25053*mg/cm3 ;       // STP
  G4Material* Nitrogen = new G4Material(name="N2", density, ncomponents=1);
  Nitrogen->AddElement(N, 2);
  
  density = 1.4289*mg/cm3 ;       // STP
  G4Material* Oxygen = new G4Material(name="O2", density, ncomponents=1);
  Oxygen->AddElement(O, 2);
  
  
  density  = 1.2928*mg/cm3 ;       // STP
  density *= 1.0e-8 ;              // pumped vacuum
  
  temperature = STP_Temperature;
  pressure = 1.0e-8*STP_Pressure;

  G4Material* LaboratoryVacuum = new G4Material(name="LaboratoryVacuum",
						density,ncomponents=3,
						kStateGas,temperature,pressure);
  LaboratoryVacuum->AddMaterial( Nitrogen, fractionmass = 0.7557 ) ;
  LaboratoryVacuum->AddMaterial( Oxygen,   fractionmass = 0.2315 ) ;
  LaboratoryVacuum->AddMaterial( Argon,    fractionmass = 0.0128 ) ;
  

  G4cout << G4endl << *(G4Material::GetMaterialTable()) << G4endl;


  // Default materials in setup.
  WorldMaterial = LaboratoryVacuum;
  GapMaterial = Fe;


  G4cout << "end material"<< G4endl;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VPhysicalVolume* PurgMagDetectorConstruction::ConstructCalorimeter()
{
  // Complete the parameters definition
  
  //The World
  WorldSizeXY  = 300.*cm;  // Cube
  WorldSizeZ   = 300.*cm;
  
  //Measurement volume
  MeasureVolumeSizeXY = 280.*cm;  // Cubic slice
  MeasureVolumeSizeZ  = 1.*cm; 

  // Position of measurement volume. 
  // SSD is Source to Surface Distance. Source in origo and measurements 50 cm 
  // below in the z-direction (symbolizin a patient at SSD = 50 cm)
 
  SSD = 50.*cm;
  MeasureVolumePosition = -(SSD + MeasureVolumeSizeZ/2); 
  

  // Geometric definition of the gap of the purging magnet. Approximation of
  // the shape of the pole gap.    

  GapSizeY1 = 10.*cm;    // length along x at the surface positioned at -dz
  GapSizeY2 = 10.*cm;    // length along x at the surface positioned at +dz
  GapSizeX1 = 10.*cm;    // length along y at the surface positioned at -dz
  GapSizeX2 = 18.37*cm;  // length along y at the surface positioned at +dz
  GapSizeZ  = 11.5*cm;   // length along z axis

  Gap1PosY = 0.*cm;
  Gap1PosX = -9.55*cm;
  Gap1PosZ = -6.89*cm;

  Gap2PosY = 0.*cm;
  Gap2PosX = 9.55*cm;
  Gap2PosZ = -6.89*cm;


  // Coordinate correction for field grif. 
  // Gap opening at z = -11.4 mm.
  // In field grid coordonates gap at z = -0.007m in field from z = 0.0m to 
  // z = 0.087m.
  // -> zOffset = -11.4-(-7) = 4.4 mm

  zOffset = 4.4*mm;  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Some out prints of the setup. 
  
  G4cout << "\n-----------------------------------------------------------"
	 << "\n      Geometry and materials"
	 << "\n-----------------------------------------------------------"
	 << "\n ---> World:" 
	 << "\n ---> " << WorldMaterial->GetName() << " in World"
	 << "\n ---> " << "WorldSizeXY: " << G4BestUnit(WorldSizeXY,"Length")
	 << "\n ---> " << "WorldSizeZ: " << G4BestUnit(WorldSizeZ,"Length");
  
#if GAP
  G4cout << "\n-----------------------------------------------------------"
	 << "\n ---> Purging Magnet:" 
	 << "\n ---> " << "Gap made of "<< GapMaterial->GetName() 
	 << "\n ---> " << "GapSizeY1: " << G4BestUnit(GapSizeY1,"Length") 
	 << "\n ---> " << "GapSizeY2: " << G4BestUnit(GapSizeY2,"Length") 
	 << "\n ---> " << "GapSizeX1: " << G4BestUnit(GapSizeX1,"Length") 
	 << "\n ---> " << "GapSizeX2: " << G4BestUnit(GapSizeX2,"Length");
#endif
  
#if MEASUREVOL
  G4cout << "\n-----------------------------------------------------------"
	 << "\n ---> Measurement Volume:" 
	 << "\n ---> " << WorldMaterial->GetName() << " in Measurement volume"
	 << "\n ---> " << "MeasureVolumeXY: " << G4BestUnit(MeasureVolumeSizeXY,"Length") 
	 << "\n ---> " << "MeasureVolumeZ: " << G4BestUnit(MeasureVolumeSizeZ,"Length")
	 << "\n ---> " << "At SSD =  " << G4BestUnit(MeasureVolumePosition,"Length");
#endif
  
  G4cout << "\n-----------------------------------------------------------\n";
  
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  //     
  // World
  //
  

  solidWorld = new G4Box("World",				       //its name
			   WorldSizeXY/2,WorldSizeXY/2,WorldSizeZ/2);  //its size
  

  logicWorld = new G4LogicalVolume(solidWorld,	        //its solid
				   WorldMaterial,	//its material
				   "World");		//its name
  
  physiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 "World",		//its name
                                 logicWorld,		//its logical volume
                                 NULL,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number

  // Visualization attributes
  G4VisAttributes* simpleWorldVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0)); //White
  simpleWorldVisAtt->SetVisibility(true);
  logicWorld->SetVisAttributes(simpleWorldVisAtt);
 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  //     
  // Measurement Volume
  //
  
#if MEASUREVOL

  solidMeasureVolume = new G4Box("MeasureVolume",				       //its name
                   MeasureVolumeSizeXY/2,MeasureVolumeSizeXY/2,MeasureVolumeSizeZ/2);  //its size

  logicMeasureVolume = new G4LogicalVolume(solidMeasureVolume,	//its solid
                                   WorldMaterial,	        //its material
                                   "MeasureVolume");		//its name
                                   
  physiMeasureVolume = new G4PVPlacement(0,			             //no rotation
  				 G4ThreeVector(0.,0.,MeasureVolumePosition), //at (0,0,0)
                                 "MeasureVolume",		             //its name
                                 logicMeasureVolume,		             //its logical volume
                                 physiWorld,			             //its mother  volume
                                 false,			                     //no boolean operation
                                 0);			                     //copy number

  // Visualization attributes
  G4VisAttributes* simpleMeasureVolumeVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0)); //White
  simpleMeasureVolumeVisAtt->SetVisibility(true);
  simpleMeasureVolumeVisAtt->SetForceSolid(true);
  logicMeasureVolume->SetVisAttributes(simpleMeasureVolumeVisAtt);

#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  //                              
  //Gap cone. Opening 20 deg. Two separate trapezoids. Iron.
  // 

#if GAP

  //Gap part 1, placed in negative x-direction.

  solidGap1 = new G4Trd("Gap1",
			GapSizeX1/2,  // Half-length along x at the surface positioned at -dz
			GapSizeX2/2,  // Half-length along x at the surface positioned at +dz
			GapSizeY1/2,  // Half-length along y at the surface positioned at -dz
			GapSizeY2/2,  // Half-length along y at the surface positioned at +dz
			GapSizeZ/2 ); // Half-length along z axis
  
  logicGap1 = new G4LogicalVolume(solidGap1,   	        //its solid
				  GapMaterial,          //its material
				  "Gap1");              //its name
  
  physiGap1 = new G4PVPlacement(0,			                    //90 deg rotation
				G4ThreeVector(Gap1PosX,Gap1PosY,Gap1PosZ),  //position
				"Gap1",		             		    //its name
				logicGap1,		                    //its logical volume
				physiWorld,		                    //its mother  volume
				false,			                    //no boolean operation
				0);			                    //copy number
  
  //Gap part 2, placed in positive x-direction.

  solidGap2 = new G4Trd("Gap2",
                	GapSizeX1/2,  // Half-length along x at the surface positioned at -dz
		        GapSizeX2/2,  // Half-length along x at the surface positioned at +dz
                	GapSizeY1/2,  // Half-length along y at the surface positioned at -dz
			GapSizeY2/2,  // Half-length along y at the surface positioned at +dz
                	GapSizeZ/2 ); // Half-length along z axis
  
  logicGap2 = new G4LogicalVolume(solidGap2,	        //its solid
				  GapMaterial,          //its material
				  "Gap2");              //its name
  
  physiGap2 = new G4PVPlacement(0,			                    //no rotation
				G4ThreeVector(Gap2PosX,Gap2PosY,Gap2PosZ),  //position
				"Gap2",		                 	    //its name
				logicGap2,		             	    //its logical volume
				physiWorld,		             	    //its mother  volume
				false,			             	    //no boolean operation
				0);			             	    //copy number

  // Visualization attributes
  G4VisAttributes* simpleGap1VisAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0)); //yellow
  simpleGap1VisAtt->SetVisibility(true);
  simpleGap1VisAtt->SetForceSolid(true);
  logicGap1->SetVisAttributes(simpleGap1VisAtt);
  
  G4VisAttributes* simpleGap2VisAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0)); //yellow
  simpleGap2VisAtt->SetVisibility(true);
  simpleGap2VisAtt->SetForceSolid(true);  
  logicGap2->SetVisAttributes(simpleGap2VisAtt);

#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  return physiWorld;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PurgMagDetectorConstruction::ConstructSDandField()
{
//  Magnetic Field - Purging magnet
//
#if MAG
  
  if (fField.Get() == 0)
    {
      //Field grid in A9.TABLE. File must be in accessible from run urn directory. 
      G4MagneticField* PurgMagField= new PurgMagTabulatedField3D("PurgMag3D.TABLE", zOffset);
      fField.Put(PurgMagField);
      
      //This is thread-local
      G4FieldManager* pFieldMgr = 
	G4TransportationManager::GetTransportationManager()->GetFieldManager();
           
      G4cout<< "DeltaStep "<<pFieldMgr->GetDeltaOneStep()/mm <<"mm" <<G4endl;
      //G4ChordFinder *pChordFinder = new G4ChordFinder(PurgMagField);

      pFieldMgr->SetDetectorField(fField.Get());
      pFieldMgr->CreateChordFinder(fField.Get());
      
    }
#endif
}
