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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:fiberMat(0),lvol_fiber(0), absorberMat(0),lvol_layer(0),
 moduleMat(0),lvol_module(0), calorimeterMat(0),lvol_calorimeter(0),
 worldMat(0),pvol_world(0), defaultMat(0)
{
  // materials
  DefineMaterials();
  
  // default parameter values of calorimeter
  //
  fiberDiameter       = 1.13*mm; 	//1.08*mm
  nbOfFibers          = 490;		//490
  distanceInterFibers = 1.35*mm;	//1.35*mm
  layerThickness      = 1.73*mm;	//1.68*mm  
  milledLayer         = 1.00*mm;    //1.40*mm ?
  nbOfLayers          = 10;		    //10
  nbOfModules         = 9;		    //9
     
  fiberLength         = (nbOfFibers+0.5)*distanceInterFibers;	//662.175*mm
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructCalorimeter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // define Elements
  //
  G4Element* H  = new G4Element("Hydrogen","H", 1,  1.01*g/mole);
  G4Element* C  = new G4Element("Carbon",  "C", 6, 12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen","N", 7, 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",  "O", 8, 16.00*g/mole);

  G4int natoms, ncomponents;
  G4double density, massfraction;				     

  // Lead
  //
  G4Material* Pb =   
  new G4Material("Lead", 82., 207.20*g/mole, density= 0.98*11.20*g/cm3);

  // Scintillator
  //
  G4Material* Sci = 
  new G4Material("Scintillator", density= 1.032*g/cm3, ncomponents=2);
  Sci->AddElement(C, natoms=8);
  Sci->AddElement(H, natoms=8);
  
  Sci->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  // Air
  //
  G4Material* Air = 
  new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, massfraction=70*perCent);
  Air->AddElement(O, massfraction=30.*perCent);

  // example of vacuum
  //
  density     = universe_mean_density;    //from PhysicalConstants.h
  G4double pressure    = 3.e-18*pascal;
  G4double temperature = 2.73*kelvin;
  G4Material* Vacuum =   
  new G4Material("Galactic", 1., 1.008*g/mole, density,
                             kStateGas,temperature,pressure);

  //attribute materials
  //
  defaultMat     = Vacuum;  
  fiberMat       = Sci;
  absorberMat    = Pb;
  moduleMat      = defaultMat;
  calorimeterMat = defaultMat;
  worldMat       = defaultMat;

  // print table
  //      
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructCalorimeter()
{
  // Cleanup old geometry
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  // fibers
  //
  G4Tubs*
  svol_fiber = new G4Tubs("fiber",			//name
                         0*mm, 0.5*fiberDiameter,	//r1, r2
			 0.5*fiberLength,		//half-length 
			 0., twopi);			//theta1, theta2
			 
  lvol_fiber = new G4LogicalVolume(svol_fiber,		//solid
                                   fiberMat,		//material
                                   "fiber");		//name
				   
  // layer
  //
  G4double sizeX = layerThickness;
  G4double sizeY = distanceInterFibers*nbOfFibers;
  G4double sizeZ = fiberLength;
  
  G4Box*      
  svol_layer = new G4Box("layer",			//name
                  0.5*sizeX, 0.5*sizeY, 0.5*sizeZ);	//size


  lvol_layer = new G4LogicalVolume(svol_layer,		//solid
                                   absorberMat,		//material
                                   "layer");		//name

  // put fibers within layer
  //
  G4double Xcenter = 0.;
  G4double Ycenter = -0.5*(sizeY + distanceInterFibers);
  
  for (G4int k=0; k<nbOfFibers; k++) {
    Ycenter += distanceInterFibers;
    new G4PVPlacement(0,		   		//no rotation
      		  G4ThreeVector(Xcenter,Ycenter,0.),    //position
                      lvol_fiber,     		   	//logical volume	
                      "fiber",	   			//name
                      lvol_layer,        		//mother
                      false,             		//no boulean operat
                      k+1);               		//copy number

  }
				   
  // modules
  //
  moduleThickness = layerThickness*nbOfLayers + milledLayer;       
  sizeX = moduleThickness;
  sizeY = fiberLength;
  sizeZ = fiberLength;
  
  G4Box*      
  svol_module = new G4Box("module",			//name
                  0.5*sizeX, 0.5*sizeY, 0.5*sizeZ);	//size

  lvol_module = new G4LogicalVolume(svol_module,	//solid
                                   absorberMat,		//material
                                   "module");		//name

  // put layers within module
  //
  Xcenter = -0.5*(nbOfLayers+1)*layerThickness;
  Ycenter =  0.25*distanceInterFibers;
  
  for (G4int k=0; k<nbOfLayers; k++) {
    Xcenter += layerThickness;
    Ycenter  = - Ycenter;
    new G4PVPlacement(0,		   		//no rotation
      		  G4ThreeVector(Xcenter,Ycenter,0.),    //position
                      lvol_layer,     		   	//logical volume	
                      "layer",	   			//name
                      lvol_module,        		//mother
                      false,             		//no boulean operat
                      k+1);               		//copy number

  }				   				   

  // calorimeter
  //
  calorThickness = moduleThickness*nbOfModules;
  sizeX = calorThickness;
  sizeY = fiberLength;
  sizeZ = fiberLength;
  
  G4Box*      
  svol_calorimeter = new G4Box("calorimeter",		//name
                  0.5*sizeX, 0.5*sizeY, 0.5*sizeZ);	//size


  lvol_calorimeter = new G4LogicalVolume(svol_calorimeter,	//solid
                                   calorimeterMat,		//material
                                   "calorimeter");		//name  

  // put modules inside calorimeter
  //  
  Xcenter = -0.5*(calorThickness + moduleThickness);
  

  for (G4int k=0; k<nbOfModules; k++) {
    Xcenter += moduleThickness;		  
    G4RotationMatrix rotm;                    //rotation matrix to place modules    
    if ((k+1)%2 == 0) rotm.rotateX(90*deg);
	G4Transform3D transform(rotm, G4ThreeVector(Xcenter,0.,0.));    
    new G4PVPlacement(transform,		   		//rotation+position
                      lvol_module,	     		//logical volume	
                      "module", 	   		    //name
                      lvol_calorimeter,        	//mother
                      false,             		//no boulean operat
                      k+1);               		//copy number
  }

  // world
  //
  sizeX = 1.2*calorThickness;
  sizeY = 1.2*fiberLength;
  sizeZ = 1.2*fiberLength;
  
  worldSizeX = sizeX;
  
  G4Box*      
  svol_world = new G4Box("world",			//name
                  0.5*sizeX, 0.5*sizeY, 0.5*sizeZ);	//size

  lvol_world = new G4LogicalVolume(svol_world,		//solid
                                   worldMat,		//material
                                   "world");		//name 
				    
  pvol_world = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 lvol_world,		//logical volume
                                 "world",		//name
                                 0,			//mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number

  //put calorimeter in world
  //  
  new G4PVPlacement(0,				//no rotation
                    G4ThreeVector(),		//at (0,0,0)
                    lvol_calorimeter,		//logical volume
                    "calorimeter",		//name
                    lvol_world,			//mother  volume
                    false,			//no boolean operation
                    0);				//copy number
		    				 
  PrintCalorParameters();
  
  // Visualization attributes
  //
  lvol_fiber->SetVisAttributes (G4VisAttributes::GetInvisible());  
  lvol_layer->SetVisAttributes (G4VisAttributes::GetInvisible());
  lvol_world->SetVisAttributes (G4VisAttributes::GetInvisible());
    
  //always return the physical World
  //
  return pvol_world;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4UnitsTable.hh"

void DetectorConstruction::PrintCalorParameters()
{
  G4cout << "\n-------------------------------------------------------------"
     << "\n ---> The calorimeter is " << nbOfModules << " Modules"
     << "\n ---> A Module is " << nbOfLayers << " Layers + 1 milled Layer";
     
  G4cout  
     << "\n ---> A Layer is " << G4BestUnit(layerThickness,"Length")  
     << " thickness of " << absorberMat->GetName();    
     
  G4cout 
     << "\n ---> A Layer includes " << nbOfFibers << " fibers of " 
     << fiberMat->GetName();
     
  G4cout 
     << "\n      ---> diameter : " << G4BestUnit(fiberDiameter,"Length")
     << "\n      ---> length   : " << G4BestUnit(fiberLength,"Length")
     << "\n      ---> distance : " << G4BestUnit(distanceInterFibers,"Length");
     
  G4cout  
     << "\n ---> The milled Layer is " << G4BestUnit(milledLayer,"Length")  
     << " thickness of " << absorberMat->GetName();
     
  G4cout 
   << "\n\n ---> Module thickness " << G4BestUnit(moduleThickness,"Length");
  
  G4cout 
   << "\n\n ---> Total calor thickness " << G4BestUnit(calorThickness,"Length")
   <<   "\n      Tranverse size        " << G4BestUnit(fiberLength,"Length");

  G4cout << "\n-------------------------------------------------------------\n";
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

void DetectorConstruction::ConstructSDandField()
{
    if ( fFieldMessenger.Get() == 0 ) {
        // Create global magnetic field messenger.
        // Uniform magnetic field is then created automatically if
        // the field value is not zero.
        G4ThreeVector fieldValue = G4ThreeVector();
        G4GlobalMagFieldMessenger* msg =
        new G4GlobalMagFieldMessenger(fieldValue);
        //msg->SetVerboseLevel(1);
        G4AutoDelete::Register(msg);
        fFieldMessenger.Put( msg );
        
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
