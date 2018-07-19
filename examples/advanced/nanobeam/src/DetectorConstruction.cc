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
// Please cite the following paper if you use this software
// Nucl.Instrum.Meth.B260:20-27, 2007
//
// Based on purging magnet advanced example.
//

#include "DetectorConstruction.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4RunManager.hh" 

// Field
#include "G4Mag_UsualEqRhs.hh"
#include "G4TransportationManager.hh"
#include "G4ClassicalRK4.hh"
#include "G4PropagatorInField.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreadLocal TabulatedField3D* DetectorConstruction::fField = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::DetectorConstruction()
{ 
 fDetectorMessenger = new DetectorMessenger(this);
 
 // Default values (square field, coef calculation, profile)
 
 fModel=1;
 fG1=-11.964623; 
 fG2=16.494652; 
 fG3=9.866770; 
 fG4=-6.244493; 
 fCoef=0; 
 fProfile=1; 
 fGrid=0;

}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

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
  
  G4double z, a;

  // Vacuum standard definition...
  density = universe_mean_density;
  G4Material* vacuum = new G4Material(name="Vacuum", z=1., a=1.01*g/mole,
	density);

  // NIST
  G4NistManager *man=G4NistManager::Instance();
  man->SetVerbose(1);

  //
  
  G4cout << G4endl << *(G4Material::GetMaterialTable()) << G4endl;

  // Default materials in setup.
  fDefaultMaterial = vacuum;
  fGridMaterial = man->FindOrBuildMaterial("G4_Ni"); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{

  fSolidWorld = new G4Box("World",		   	//its name
			   12*m/2,12*m/2,22*m/2);  	//its size
  

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


  // MAGNET VOLUME 

  fSolidVol = new G4Box("Vol",				//its name
			   10*m/2,10*m/2,9.120*m/2);  	//its size
  

  fLogicVol = new G4LogicalVolume(fSolidVol,	        //its solid
				  fDefaultMaterial,	//its material
				  "Vol");		//its name
  
  fPhysiVol = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(0,0,-4310*mm),	//at (0,0,0)
                                 "Vol",			//its name
                                 fLogicVol,		//its logical volume
                                 fPhysiWorld,		//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number

  // GRID
  
  if (fGrid==1)
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

  fSolidGridVol= new G4Box("GridVolume",x_grid,y_grid,z_grid);   //its size
  
  fLogicGridVol = new G4LogicalVolume(fSolidGridVol,  		//its solid
				      fGridMaterial,            //its material
				      "GridVolume");		//its name
  
  fPhysiGridVol = new G4PVPlacement(0,				//no rotation
  				 G4ThreeVector(0,0,grid_Zpos),	// origin
                                 fLogicGridVol,			//its logical volume
                                 "GridVolume",			//its name
                                 fLogicWorld,	        	//its mother  volume
                                 false,				//no boolean operation
                                 0);	

  // Holes in grid
  
  G4double holeSize= 9e-3*mm;
  G4double pix_grid=1.3e-2*mm;
  G4int    num_half_grid=100;

  fSolidGridVol_Hole= new G4Box("GridHole",holeSize/2,holeSize/2,z_grid);   //its size
  
  fLogicGridVol_Hole = new G4LogicalVolume(fSolidGridVol_Hole,  	    //its solid
				   fDefaultMaterial,                        //its material
				   "GridHole");		                    //its name

 
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

   		fPhysiGridVol_Hole  = new G4PVPlacement(0,		//no rotation
			  	 G4ThreeVector(x0_grid,y0_grid,z0_grid),//origin
                                 fLogicGridVol_Hole,			//its logical volume
  			         "GridHole",				//its name
                                 fLogicGridVol,	        		//its mother  volume
                                 false,					//no boolean operation
                                 number_index_grid);
	}	
  }

  // Grid imaging plane
  
  G4double ContVolSizeXY = 1*m;
  G4double ImPlaneWidth = 0.001*mm;
 
  fSolidControlVol_GridShadow =
    new G4Box
    ("ControlVol_GridShadow", ContVolSizeXY/2, ContVolSizeXY/2 , ImPlaneWidth/2);
 
  fLogicControlVol_GridShadow = 
    new G4LogicalVolume
    (fSolidControlVol_GridShadow, fDefaultMaterial, "ControlVol_GridShadow");
  
  fPhysiControlVol_GridShadow = 
    new G4PVPlacement 
    ( 0, G4ThreeVector(0,0,(250+300)*mm), fLogicControlVol_GridShadow, "ControlVol_GridShadow",
      fLogicWorld, false, 0);
     
 
  } // end GRID
  
  // STEP MINIMUM SIZE 
  fLogicVol->SetUserLimits(new G4UserLimits(1*mm));

  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::SetG1(G4float value)
{
  fG1 = value; 
  G4RunManager::GetRunManager()->ReinitializeGeometry();   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::SetG2(G4float value)
{
  fG2 = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::SetG3(G4float value)
{
  fG3 = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::SetG4(G4float value)
{
  fG4 = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::SetModel(G4int modelChoice)
{
  if (modelChoice==1) fModel=1;
  if (modelChoice==2) fModel=2;
  if (modelChoice==3) fModel=3;
  G4RunManager::GetRunManager()->ReinitializeGeometry();   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::SetCoef(G4int val)
{
  fCoef=val;
  G4RunManager::GetRunManager()->ReinitializeGeometry();   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int DetectorConstruction::GetCoef()
{
  return fCoef;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::SetProfile(G4int myProfile)
{
  fProfile=myProfile;
  G4RunManager::GetRunManager()->ReinitializeGeometry();   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::SetGrid(G4int myGrid)
{
  fGrid=myGrid;
  G4RunManager::GetRunManager()->ReinitializeGeometry();   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::ConstructSDandField()
{
      fField = new TabulatedField3D(fG1, fG2, fG3, fG4, fModel); 
   
      //This is thread-local
      G4FieldManager* fFieldMgr = 
	G4TransportationManager::GetTransportationManager()->GetFieldManager();
           
      G4Mag_UsualEqRhs* fEquation = new G4Mag_UsualEqRhs (fField);

      G4ClassicalRK4* fStepper = new G4ClassicalRK4 (fEquation);

      G4ChordFinder* fChordFinder = new G4ChordFinder(fField,1e-9*m,fStepper);

      fFieldMgr->SetChordFinder(fChordFinder);
      fFieldMgr->SetDetectorField(fField);    
 
      fFieldMgr->GetChordFinder()->SetDeltaChord(1e-9*m);
      fFieldMgr->SetDeltaIntersection(1e-9*m);
      fFieldMgr->SetDeltaOneStep(1e-9*m);     
      
      // To avoid G4MagIntegratorDriver::OneGoodStep:Stepsize underflows in Stepper
      
      if (fCoef==1)
      {
        G4PropagatorInField* fPropInField =
          G4TransportationManager::GetTransportationManager()->GetPropagatorInField();
        fPropInField->SetMinimumEpsilonStep(1e-11);
        fPropInField->SetMaximumEpsilonStep(1e-10); 

      } 
      else
      {
        G4PropagatorInField* fPropInField =
          G4TransportationManager::GetTransportationManager()->GetPropagatorInField();
        fPropInField->SetMinimumEpsilonStep(1e-9);
        fPropInField->SetMaximumEpsilonStep(1e-8);
      }

}

