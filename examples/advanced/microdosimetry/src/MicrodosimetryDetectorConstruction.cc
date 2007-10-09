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
// $Id: MicrodosimetryDetectorConstruction.cc,v 1.1 2007-10-09 08:00:28 sincerti Exp $
// -------------------------------------------------------------------

#include "MicrodosimetryDetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrodosimetryDetectorConstruction::MicrodosimetryDetectorConstruction()
  
  
{
  WorldSizeXY=WorldSizeZ=0;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrodosimetryDetectorConstruction::~MicrodosimetryDetectorConstruction()
{
  delete defaultMaterial;
  delete KgmMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* MicrodosimetryDetectorConstruction::Construct()
  
{
  DefineMaterials();
  return ConstructMicrodosimetryLine();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicrodosimetryDetectorConstruction::DefineMaterials()
{ 

  G4String name, symbol;             
  G4double density;            
  
  G4int ncomponents, natoms;
  G4double z, a;
  
  // Define Elements 
  G4Element*   H  = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);
  G4Element*   O  = new G4Element ("Oxygen"  , "O", 8. , 16.00*g/mole);
 
  // Vaccum standard definition...
  density = universe_mean_density;
  G4Material* vacuum = new G4Material(name="Vacuum", z=1., a=1.01*g/mole,
				      density);  
  // Water 
  density = 1.000*g/cm3;
  G4Material* H2O = new G4Material(name="H2O"  , density, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
 
 
  // Materials in setup.
  defaultMaterial 	= vacuum;
  KgmMaterial 		= H2O;
  
  // DISPLAY MATERIALS
  G4cout << G4endl << *(G4Material::GetMaterialTable()) << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* MicrodosimetryDetectorConstruction::ConstructMicrodosimetryLine()
{
  // WORLD
  WorldSizeXY  = 5*cm;
  WorldSizeZ   = 6*mm;
   
  //*************
  // WORLD VOLUME
  //*************
  
  solidWorld = new G4Box("World",				       //its name
			 WorldSizeXY/2,WorldSizeXY/2,WorldSizeZ/2);  //its size
  
  
  logicWorld = new G4LogicalVolume(solidWorld,	        //its solid
				   KgmMaterial,	//its material
				   "World");		//its name
  
  physiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 "World",		//its name
                                 logicWorld,		//its logical volume
                                 NULL,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number

  //****
  // KGM   
  //****
    
  solidKgm = new G4Box("KGM", 2.*cm, 2.*cm, 3*mm/2);

  logicKgm = new G4LogicalVolume(solidKgm, KgmMaterial, "KGM");
  
  physiKgm = new G4PVPlacement(0,
			       G4ThreeVector(0,0,1.5*mm),
			       "KGM",
			       logicKgm,
			       physiWorld, 
			       false,
			       0);

  // ************
  // CELL PHANTOM
  // ************

  solidPhantom = new G4Box("Phantom", 
  	myMicrodosimetryPhantomConfiguration.GetPixelSizeX()/2, 
	myMicrodosimetryPhantomConfiguration.GetPixelSizeY()/2, 
	myMicrodosimetryPhantomConfiguration.GetPixelSizeZ()/2); 
  
  logicPhantom = new G4LogicalVolume(solidPhantom,KgmMaterial,"Phantom",0,0,0);
    
    // PHANTOM MASSES

  SetNbOfPixelsInPhantom (myMicrodosimetryPhantomConfiguration.GetPhantomTotalPixels());

  SetMassNucleus(myMicrodosimetryPhantomConfiguration.GetNucleusMass());

  SetMassCytoplasm(myMicrodosimetryPhantomConfiguration.GetCytoplasmMass());

  // PHANTOM

  phantomParam = new MicrodosimetryCellParameterisation
  	(myMicrodosimetryPhantomConfiguration.GetPhantomTotalPixels(),
	 myMicrodosimetryPhantomConfiguration.GetPixelSizeX()/2,
	 myMicrodosimetryPhantomConfiguration.GetPixelSizeY()/2,
	 myMicrodosimetryPhantomConfiguration.GetPixelSizeZ()/2,
	 KgmMaterial,KgmMaterial,
	 KgmMaterial,KgmMaterial,
	 KgmMaterial,KgmMaterial
	 );

  physiPhantom = new G4PVParameterised(
                            "Phantom",       // their name
                            logicPhantom,    // their logical volumr
//                            logicCyto,       // Mother logical volume
                            logicKgm,       // Mother logical volume
			    kZAxis,          // Are placed along this axis 
                            phantomParam->GetNoBoxes(),    // Number of boxes
                            phantomParam,false);   // The parametrisation

  G4cout << " ==========> The phantom contains " 
    << myMicrodosimetryPhantomConfiguration.GetPhantomTotalPixels() << " voxels " << G4endl;		    
  G4cout << G4endl; 
				    		    
  // USER LIMITS ON STEP LENGTH
  
  // VISUALISATION ATTRIBUTES (for phantom, see in Parameterisation class)
  
  G4VisAttributes* simpleBoxAttKGM= new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  simpleBoxAttKGM->SetDaughtersInvisible(false);
  simpleBoxAttKGM->SetForceSolid(false);
  
  
  logicKgm->SetVisAttributes(simpleBoxAttKGM);

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
  return physiWorld;
}
