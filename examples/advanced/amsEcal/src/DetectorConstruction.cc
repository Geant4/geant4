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
// $Id: DetectorConstruction.cc,v 1.1 2009-04-16 11:05:40 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4RotationMatrix.hh"

#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UnitsTable.hh"

#include "G4VisAttributes.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:fiberMat(0),lvol_fiber(0), absorberMat(0),lvol_layer(0),
 superLayerMat(0),lvol_superlayer(0), calorimeterMat(0),lvol_calorimeter(0),
 worldMat(0),pvol_world(0), defaultMat(0), magField(0)
{
  // materials
  DefineMaterials();
  
  // default parameter values of calorimeter
  //
  fiberDiameter       = 1.08*mm;	//1.08*mm
  nbOfFibers          = 490;		//490
  distanceInterFibers = 1.35*mm;	//1.35*mm
  distanceInterLayers = 1.68*mm;	//1.68*mm
  nbOfLayers          = 10;		//10
  nbOfSuperLayers     = 9;		//9 
  
  nbOfLayersPerPixel  = 5;              //5 (must divide nbOfLayers)
  
  fiberLength         = (nbOfFibers+1)*distanceInterFibers;	//658*mm    
  nbOfPixels          = (nbOfLayers*nbOfSuperLayers)/nbOfLayersPerPixel;  //18
  
  // create commands for interactive definition of the calorimeter
  detectorMessenger = new DetectorMessenger(this);   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
  delete detectorMessenger;
}

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
  new G4Material("Lead", 82., 207.20*g/mole, density= 11.20*g/cm3);

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
  superLayerMat  = defaultMat;
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
  G4double layerThick = distanceInterLayers;
  G4double sizeX = layerThick;
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
                      k);               		//copy number

  }
				   
  // super layer
  //
  G4double superLayerThick = layerThick*nbOfLayers;
  sizeX = superLayerThick;
  sizeY = fiberLength;
  sizeZ = fiberLength;
  
  G4Box*      
  svol_superlayer = new G4Box("superLayer",		//name
                  0.5*sizeX, 0.5*sizeY, 0.5*sizeZ);	//size

  lvol_superlayer = new G4LogicalVolume(svol_superlayer,	//solid
                                   absorberMat,			//material
                                   "superLayer");		//name

  // put layers within superlayer
  //
  Xcenter = -0.5*(superLayerThick + layerThick);
  Ycenter =  0.25*distanceInterFibers;
  
  for (G4int k=0; k<nbOfLayers; k++) {
    Xcenter += layerThick;
    Ycenter  = - Ycenter;
    new G4PVPlacement(0,		   		//no rotation
      		  G4ThreeVector(Xcenter,Ycenter,0.),    //position
                      lvol_layer,     		   	//logical volume	
                      "layer",	   			//name
                      lvol_superlayer,        		//mother
                      false,             		//no boulean operat
                      k);               		//copy number

  }				   				   

  // calorimeter
  //
  G4double calorimeterThick = superLayerThick*nbOfSuperLayers;
  sizeX = calorimeterThick;
  sizeY = fiberLength;
  sizeZ = fiberLength;
  
  G4Box*      
  svol_calorimeter = new G4Box("calorimeter",		//name
                  0.5*sizeX, 0.5*sizeY, 0.5*sizeZ);	//size


  lvol_calorimeter = new G4LogicalVolume(svol_calorimeter,	//solid
                                   calorimeterMat,		//material
                                   "calorimeter");		//name  

  // put superLayers inside calorimeter
  //  
  Xcenter = -0.5*(calorimeterThick + superLayerThick);
  
  //rotation matrix to place superLayers
  G4RotationMatrix* rotm = 0;  
  G4RotationMatrix* rotmX = new G4RotationMatrix();
  rotmX->rotateX(90*deg);
    
  for (G4int k=0; k<nbOfSuperLayers; k++) {
    rotm = 0;
    if ((k+1)%2 == 0) rotm = rotmX;
    Xcenter += superLayerThick;    
    new G4PVPlacement(rotm,		   		//rotation
      		  G4ThreeVector(Xcenter,0.,0.),		//position
                      lvol_superlayer,     		//logical volume	
                      "superLayer",	   		//name
                      lvol_calorimeter,        		//mother
                      false,             		//no boulean operat
                      k);               		//copy number
  }

  // world
  //
  sizeX = 1.2*calorimeterThick;
  sizeY = 1.2*fiberLength;
  sizeZ = 1.2*fiberLength;
  
  worldSizeX = sizeX;
  
  G4Box*      
  svol_world = new G4Box("world",			//name
                  0.5*sizeX, 0.5*sizeY, 0.5*sizeZ);	//size

  G4LogicalVolume*
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
  lvol_fiber->SetVisAttributes (G4VisAttributes::Invisible);  
  lvol_layer->SetVisAttributes (G4VisAttributes::Invisible);
  lvol_world->SetVisAttributes (G4VisAttributes::Invisible);
  
  //always return the physical World
  //
  return pvol_world;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintCalorParameters()
{
  G4cout << "\n-------------------------------------------------------------"
     << "\n ---> The calorimeter is " << nbOfSuperLayers << " superLayers"
     << "\n ---> A superLayer is " << nbOfLayers << " Layers";
     
  G4cout 
     << "\n ---> A Layer is " << G4BestUnit(distanceInterLayers,"Length")  
     << " thickness of " << absorberMat->GetName();    
     
  G4cout 
     << "\n ---> A Layer include " << nbOfFibers << " fibers of " 
     << fiberMat->GetName();
     
  G4cout 
     << "\n      ---> diameter : " << G4BestUnit(fiberDiameter,"Length")
     << "\n      ---> length   : " << G4BestUnit(fiberLength,"Length")
     << "\n      ---> distance : " << G4BestUnit(distanceInterFibers,"Length");

  G4cout << "\n-------------------------------------------------------------\n";
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

void DetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along Z axis
  //
  G4FieldManager* fieldMgr
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  if(magField) delete magField;		//delete the existing magn field

  if(fieldValue!=0.)			// create a new one if non nul
  { magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);
  } else {
    magField = 0;
    fieldMgr->SetDetectorField(magField);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructCalorimeter());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
