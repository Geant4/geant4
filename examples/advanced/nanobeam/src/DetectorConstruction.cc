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
// $Id: DetectorConstruction.cc,v 1.3 2008-12-18 12:56:24 gunter Exp $
// -------------------------------------------------------------------

#include "DetectorConstruction.hh"
#include "G4NistManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::DetectorConstruction()
{ 
 detectorMessenger = new DetectorMessenger(this);
 gradientsInitialized=false;
 G1=0; G2=0; G3=0; G4=0; coef=0; profile=0; grid=0;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::~DetectorConstruction()
{ delete detectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* DetectorConstruction::Construct()

{
  DefineMaterials();
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::DefineMaterials()
{ 
  G4String name, symbol;             
  G4double density;            
  
  G4int ncomponents, natoms;
  G4double z, a;
  
  // Define Elements  
  G4Element*   H  = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);
  G4Element*   O  = new G4Element ("Oxygen"  , "O", 8. , 16.00*g/mole);

  // Water 
  density = 1.000*g/cm3;
  G4Material* H2O = new G4Material(name="H2O"  , density, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);

  // Vacuum standard definition...
  density = universe_mean_density;
  G4Material* vacuum = new G4Material(name="Vacuum", z=1., a=1.01*g/mole,
	density);

  // NIST
  G4NistManager *man=G4NistManager::Instance();
  man->SetVerbose(1);

  G4cout << G4endl << *(G4Material::GetMaterialTable()) << G4endl;

  // Default materials in setup.
  defaultMaterial = vacuum;
  waterMaterial = H2O;
  gridMaterial = man->FindOrBuildMaterial("G4_Ni"); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{

  static G4bool fieldIsInitialized = false;
  if(!fieldIsInitialized && gradientsInitialized)
  {
      G4FieldManager* pFieldMgr;
      G4MagIntegratorStepper* pStepper;
      G4Mag_UsualEqRhs* pEquation;
    
      G4MagneticField* Field= new TabulatedField3D(G1, G2, G3, G4, model);
      
      pEquation = new G4Mag_UsualEqRhs (Field);
      pStepper = new G4ClassicalRK4 (pEquation);
      pFieldMgr=G4TransportationManager::GetTransportationManager()->GetFieldManager();
      
      G4ChordFinder *pChordFinder = new G4ChordFinder(Field,1e-9*m,pStepper);
      pFieldMgr->SetChordFinder( pChordFinder );
      
      pFieldMgr->SetDetectorField(Field);
      
      fieldIsInitialized = true;
      
      // tuned parameters
      pFieldMgr->GetChordFinder()->SetDeltaChord(1.e-9*m);
      pFieldMgr->SetDeltaIntersection(1.e-9*m);
      pFieldMgr->SetDeltaOneStep(1.e-9*m);     

      G4PropagatorInField *propInField;
      propInField =
       G4TransportationManager::GetTransportationManager()->GetPropagatorInField();
      propInField->SetMinimumEpsilonStep(1e-11);
      propInField->SetMaximumEpsilonStep(1e-10);

    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  solidWorld = new G4Box("World",		   	//its name
			   12*m/2,12*m/2,22*m/2);  	//its size
  

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


  // MAGNET VOLUME 

  solidVol = new G4Box("Vol",				//its name
			   10*m/2,10*m/2,9.120*m/2);  	//its size
  

  logicVol = new G4LogicalVolume(solidVol,	        //its solid
				   defaultMaterial,	//its material
				   "Vol");		//its name
  
  physiVol = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(0,0,-4310*mm),	//at (0,0,0)
                                 "Vol",			//its name
                                 logicVol,		//its logical volume
                                 physiWorld,		//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number

  // GRID
  
  if (grid==1)
  {
  
  G4cout << G4endl;
  
  G4cout << " ********************** " << G4endl;
  G4cout << " **** GRID IN PLACE *** " << G4endl;
  G4cout << " ********************** " << G4endl;

  G4double x_grid=5.0*mm;    
  G4double y_grid=5.0*mm;
  G4double grid_Zpos=(250+200)*mm;      // 250+10 mm for object size of 50µm diam
  //G4double thickness_grid=10*micrometer;
  G4double thickness_grid=100*micrometer;
  G4double z_grid=thickness_grid/2.0; 

  solidGridVol= new G4Box("GridVolume",x_grid,y_grid,z_grid);   //its size
  
  logicGridVol = new G4LogicalVolume(solidGridVol,  		//its solid
				   gridMaterial,               	//its material
				   "GridVolume");		//its name
  
  physiGridVol = new G4PVPlacement(0,				//no rotation
  				 G4ThreeVector(0,0,grid_Zpos),	// origin
                                 logicGridVol,			//its logical volume
                                 "GridVolume",			//its name
                                 logicWorld,	        	//its mother  volume
                                 false,				//no boolean operation
                                 0);	

  // Holes in grid
  
  G4double holeSize= 9e-3*mm;
  G4double pix_grid=1.3e-2*mm;
  G4int    num_half_grid=100;

  solidGridVol_Hole= new G4Box("GridHole",holeSize/2,holeSize/2,z_grid);   //its size
  
  logicGridVol_Hole = new G4LogicalVolume(solidGridVol_Hole,  	    //its solid
				   defaultMaterial,                 //its material
				   "GridHole");		            //its name

 
  for(int i=-num_half_grid;i<num_half_grid;i++)
  {
    	for (int j=-num_half_grid;j<num_half_grid;j++)
	{

    		G4double  x0_grid,y0_grid,z0_grid;
    		G4int  number_index_grid;

    		x0_grid=pix_grid*i;
    		y0_grid=pix_grid*j;
    		z0_grid=0.0*mm;

		number_index_grid=(i+num_half_grid)*1000+(j+num_half_grid);

   		physiGridVol_Hole  = new G4PVPlacement(0,		//no rotation
			  	 G4ThreeVector(x0_grid,y0_grid,z0_grid),//origin
                                 logicGridVol_Hole,			//its logical volume
  			         "GridHole",				//its name
                                 logicGridVol,	        		//its mother  volume
                                 false,					//no boolean operation
                                 number_index_grid);
	}	
  }

  // Grid imaging plane
  
  G4double ContVolSizeXY = 1*m;
  G4double ImPlaneWidth = 0.001*mm;
 
  solidControlVol_GridShadow =
    new G4Box
    ("ControlVol_GridShadow", ContVolSizeXY/2, ContVolSizeXY/2 , ImPlaneWidth/2);
 
  logicControlVol_GridShadow = 
    new G4LogicalVolume
    (solidControlVol_GridShadow, defaultMaterial, "ControlVol_GridShadow");
  
  physiControlVol_GridShadow = 
    new G4PVPlacement
    //( 0, G4ThreeVector(0,0,(250+250)*mm), logicControlVol_GridShadow, "ControlVol_GridShadow",logicWorld, false, 0);
 
    ( 0, G4ThreeVector(0,0,(250+300)*mm), logicControlVol_GridShadow, "ControlVol_GridShadow",logicWorld, false, 0);
     
 
  } // end GRID
  
  // STEP MINIMUM SIZE 
  logicVol->SetUserLimits(new G4UserLimits(1*mm));

  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::SetG1(G4float value)
{
  G1 = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::SetG2(G4float value)
{
  G2 = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::SetG3(G4float value)
{
  G3 = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::SetG4(G4float value)
{
  G4 = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::SetModel(G4int modelChoice)
{
if (modelChoice==1) model=1;
if (modelChoice==2) model=2;
if (modelChoice==3) model=3;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4RunManager.hh" 
 
void DetectorConstruction::UpdateGeometry()
{
  gradientsInitialized=true;
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::SetCoef()
{
  coef=1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int DetectorConstruction::GetCoef()
{
  return coef;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::SetProfile(G4int myProfile)
{
  profile=myProfile;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::SetGrid(G4int myGrid)
{
  grid=myGrid;
}

