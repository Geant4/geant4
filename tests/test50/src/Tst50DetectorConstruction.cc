//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: Tst50DetectorConstruction.cc,v 1.3 2002-12-16 13:50:08 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#include "Tst50DetectorConstruction.hh"
#include "Tst50DetectorMessenger.hh"
#include "G4UserLimits.hh"

#include "Tst50TrackerSD.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SDManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "Tst50TrackerSD.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
Tst50DetectorConstruction::Tst50DetectorConstruction()
:solidWorld(0),  logicWorld(0),  physiWorld(0),
 solidTarget(0), logicTarget(0), physiTarget(0), 
 solidTracker(0),logicTracker(0),physiTracker(0), 
 TargetMater(0), ChamberMater(0)
{
  //  fpMagField = new Tst50MagneticField();
  detectorMessenger = new Tst50DetectorMessenger(this);

// set fUserLimit true to have the limit on step lenght
 
 fUseUserLimits = false;//non fa nulla se c'e false ,devo mettere true 
 theMaxStep= 0.0001*mm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
Tst50DetectorConstruction::~Tst50DetectorConstruction()
{
  // delete fpMagField;
  delete detectorMessenger;             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* Tst50DetectorConstruction::Construct()
{
//--------- Material definition ---------

  G4double a, iz, z, density;
  G4String name, symbol;
  G4double temperature, pressure;
  G4int nel;

  //Air
    a = 14.01*g/mole;
    G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a);
    a = 16.00*g/mole;
    G4Element* elO = new G4Element(name="Oxigen", symbol="O", iz=8., a);
    density = 1.29*mg/cm3;
    G4Material* Air = new G4Material(name="Air", density, nel=2);
    Air->AddElement(elN, .7);
    Air->AddElement(elO, .3);

 
  density     = universe_mean_density;
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  a=1.01*g/mole;
  G4Material*Vacuum = new G4Material("Galactic", z=1., a,
                                    density,kStateGas,temperature,pressure);


 //Cs  http://pearl1.lanl.gov/periodic/elements/55.html
    a = 132.9054*g/mole;
    G4Element* Cs = new G4Element(name="Cesium",symbol="Cs", z=55., a);
 
    //I  http://pearl1.lanl.gov/periodic/elements/53.html
    a = 126.9045*g/mole;
   G4Element* I = new G4Element(name="Iodium",symbol="I" , z=53., a);
 
 
density=4.51*g/cm3;
 G4Material* CsI=new G4Material(name="CsI", density, nel=2);

    CsI->AddElement(Cs, .5);
    CsI->AddElement(I, .5);
  
a=26.98*g/mole;
  density=2.699*g/cm3;
  G4Material* elAll=new G4Material("Aluminum", z=13.,a,density);

 
//Pb
    a = 207.19*g/mole;
    density = 11.35*g/cm3;
    G4Material* Pb = new G4Material(name="Pb", z=82., a, density);
    
  //Xenon gas
    density     = 5.458*mg/cm3;    
    pressure    = 1*atmosphere;
    temperature = 293.15*kelvin;
    G4Material* Xenon = new G4Material(name="XenonGas", z=54., a=131.29*g/mole,
                        density, kStateGas,temperature,pressure);

  // Print all the materials defined.
  //
  

//--------- Sizes of the principal geometrical components (solids)  ---------
  
    G4double xworld = 2.0*m;
    G4double yworld = 2.0*m;
    G4double zworld = 4.0*m;

    G4double xtarget = 20.0*cm;
    G4double ytarget = 20.0*cm;
    G4double ztarget =0.5*mm;
      
   

//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
  
  //------------------------------ 
  // World
  //------------------------------ 


 solidWorld= new G4Box("world", xworld, yworld, zworld);
 logicWorld= new G4LogicalVolume( solidWorld, Vacuum, "World", 0, 0, 0);
  
  //  Must place the World Physical volume unrotated at (0,0,0).
  // 
  physiWorld = new G4PVPlacement(0,               // no rotation
                                 G4ThreeVector(), // at (0,0,0)
				 "World",         // its name
                                 logicWorld,      // its logical volume
                                 0,               // its mother  volume
                                 false,           // no boolean operations
                                 0);              // no field specific to volume
				 
  //------------------------------ 
  // Target
  //------------------------------
  
  //G4ThreeVector positionTarget = G4ThreeVector(0,0,-(targetSize+trackerSize));
   
  solidTarget = new G4Box("target", xtarget, ytarget, ztarget);
  //  logicTarget = new G4LogicalVolume(solidTarget,TargetMater,"Target",0,0,0);
  logicTarget = new G4LogicalVolume(solidTarget,CsI,"Target",0,0,0);
 
 // create UserLimits
  if (theUserLimits != NULL) delete theUserLimits;
  theUserLimits= new G4UserLimits(//DBL_MAX,  //step max
					      //DBL_MAX,  // track max
					      theMaxStep);
	
 // attach UserLimits   
  if (fUseUserLimits) {
    logicTarget->SetUserLimits(theUserLimits);
  }

  physiTarget = new G4PVPlacement(0,               // no rotation
				  G4ThreeVector(),  // at (x,y,z)
				  "Target",        // its name
				  logicTarget,     // its logical volume
				  physiWorld,      // its mother  volume
				  false,           // no boolean operations
				  0);              // no particular field 

  //G4cout << "Target is " << fTargetLength/cm << " cm of " 
  //       << TargetMater->GetName() << G4endl;

  //------------------------------ 
  // Tracker
  //------------------------------
  
  //G4ThreeVector positionTracker = G4ThreeVector(0,0,0);
  
  
  //------------------------------------------------ 
  // Sensitive detectors
  //------------------------------------------------ 

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  if (pTargetSD ==0)
    {
  G4String targetSD_name = "target";
  pTargetSD = new Tst50TrackerSD( targetSD_name  );
  if(pTargetSD)
   { SDman->AddNewDetector( pTargetSD );
  logicTarget->SetSensitiveDetector( pTargetSD );}
    }
//--------- Visualization attributes -------------------------------

  G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  logicWorld  ->SetVisAttributes(BoxVisAtt);  
  logicTarget ->SetVisAttributes(BoxVisAtt);
 
  
 
  
//--------- example of User Limits -------------------------------

  // below is an example of how to set tracking constraints in a given
  // logical volume(see also in N02PhysicsList how to setup the process
  // G4UserSpecialCuts).  
  // Sets a max Step length in the tracker region
  // G4double maxStep = 0.5*ChamberWidth, maxLength = 2*fTrackerLength;
  // G4double maxTime = 0.1*ns, minEkin = 10*MeV;
  // logicTracker->SetUserLimits(new G4UserLimits(maxStep,maxLength,maxTime,
  //                                               minEkin));
  
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void Tst50DetectorConstruction::setTargetMaterial(G4String materialName)
{
  // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialName);  
  if (pttoMaterial)
     {TargetMater = pttoMaterial;
      logicTarget->SetMaterial(pttoMaterial); 
      //G4cout << "\n----> The target is " << fTargetLength/cm << " cm of "
      //       << materialName << G4endl;
     }             
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 void  Tst50DetectorConstruction::SetMaxStep(G4double value)
{ 
  theMaxStep = value; 
  if (theUserLimits != NULL) 
  {
    theUserLimits->SetMaxAllowedStep(value);
  }
}

void  Tst50DetectorConstruction::UseUserLimits(G4bool isUse) 
{
  fUseUserLimits = isUse;
  if ( fUseUserLimits && (theUserLimits!= NULL)) 
  {logicTarget->SetUserLimits(theUserLimits);
  }    
} 












