// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02DetectorConstruction.cc,v 1.6 2001-01-31 18:05:28 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..... 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
 
#include "ExN02DetectorConstruction.hh"
#include "ExN02DetectorMessenger.hh"
#include "ExN02ChamberParameterisation.hh"
#include "ExN02MagneticField.hh"
#include "ExN02TrackerSD.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SDManager.hh"

///#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
 
ExN02DetectorConstruction::ExN02DetectorConstruction()
:solidWorld(NULL),  logicWorld(NULL),  physiWorld(NULL),
 solidTarget(NULL), logicTarget(NULL), physiTarget(NULL), 
 solidTracker(NULL),logicTracker(NULL),physiTracker(NULL), 
 solidChamber(NULL),logicChamber(NULL),physiChamber(NULL), 
 TargetMater(NULL), ChamberMater(NULL),fpMagField(NULL),
 fWorldLength(0.),  fTargetLength(0.), fTrackerLength(0.),
 NbOfChambers(0) ,  ChamberWidth(0.),  ChamberSpacing(0.)
{
  fpMagField = new ExN02MagneticField();
  detectorMessenger = new ExN02DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
 
ExN02DetectorConstruction::~ExN02DetectorConstruction()
{ 
  delete fpMagField;
  delete detectorMessenger;             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
 
G4VPhysicalVolume* ExN02DetectorConstruction::Construct()
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

  //Al
    a = 26.98*g/mole;
    density = 2.7*g/cm3;
    G4Material* Aluminium = new G4Material(name="Al", z=13., a, density);

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
    G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;

//--------- Sizes of the principal geometrical components (solids)  ---------
  
  NbOfChambers = 5;
  ChamberWidth = 20*cm;
  ChamberSpacing = 80*cm;
  
  fTrackerLength = (NbOfChambers+1)*ChamberSpacing; // Full length of the Tracker
  fTargetLength  = 5.0 * cm;                        // Full length of the Target
  
  TargetMater  = Pb;
  ChamberMater = Xenon;
  
  fWorldLength= 1.2 *(fTargetLength+fTrackerLength);
   
  G4double targetSize  = 0.5*fTargetLength;    // Half length of the Target  
  G4double trackerSize = 0.5*fTrackerLength;   // Half length of the Tracker
      
//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
  
  //------------------------------ 
  // World
  //------------------------------ 

  G4double HalfWorldLength = 0.5*fWorldLength;
  
  solidWorld= new G4Box("World",HalfWorldLength,HalfWorldLength,HalfWorldLength);
  logicWorld= new G4LogicalVolume( solidWorld, Air, "World", 0, 0, 0);
  
  //  Must place the World Physical volume unrotated at (0,0,0).
  // 
  physiWorld = new G4PVPlacement(0,               // no rotation
                                 G4ThreeVector(), // at (0,0,0)
				 "WorldPV",       // its name
                                 logicWorld,      // its logical volume
                                 0,               // its mother  volume
                                 false,           // no boolean operations
                                 0);              // no field specific to volume
				 
  //------------------------------ 
  // Target
  //------------------------------
  
  G4ThreeVector positionTarget = G4ThreeVector(0,0,-(targetSize+trackerSize));
   
  solidTarget = new G4Box("TargetSolid",targetSize,targetSize,targetSize);
  logicTarget = new G4LogicalVolume(solidTarget,TargetMater,"TargetLV",0,0,0);
  physiTarget = new G4PVPlacement(0,               // no rotation
				  positionTarget,  // at (x,y,z)
				  "TargetPV",      // its name
				  logicTarget,     // its logical volume
				  physiWorld,      // its mother  volume
				  false,           // no boolean operations
				  0);              // no particular field 

  G4cout << "Target is " << fTargetLength/cm << " cm of " 
         << TargetMater->GetName() << G4endl;

  //------------------------------ 
  // Tracker
  //------------------------------
  
  G4ThreeVector positionTracker = G4ThreeVector(0,0,0);
  
  solidTracker = new G4Box("Tracker",trackerSize,trackerSize,trackerSize);
  logicTracker = new G4LogicalVolume(solidTracker , Air, "Tracker",0,0,0);  
  physiTracker = new G4PVPlacement(0,              // no rotation
				  positionTracker, // at (x,y,z)
				  "Tracker",       // its name
				  logicTracker,    // its logical volume
				  physiWorld,      // its mother  volume
				  false,           // no boolean operations
				  0);              // no particular field 

  //------------------------------ 
  // Tracker segments
  //------------------------------
  // 
  // An example of Parameterised volumes
  // dummy values for G4Box -- modified by parameterised volume

  solidChamber = new G4Box("chamber", 100*cm, 100*cm, 10*cm); 
  logicChamber = new G4LogicalVolume(solidChamber,ChamberMater,"Chamber",0,0,0);
  
  G4double firstPosition = -trackerSize + 0.5*ChamberWidth;
  G4double firstLength = fTrackerLength/10;
  G4double lastLength  = fTrackerLength;
   
  G4VPVParameterisation* chamberParam = new ExN02ChamberParameterisation(  
			   NbOfChambers,          // NoChambers 
			   firstPosition,         // Z of center of first 
			   ChamberSpacing,        // Z spacing of centers
			   ChamberWidth,          // Width Chamber 
			   firstLength,           // lengthInitial 
			   lastLength);           // lengthFinal
			   
  // dummy value : kZAxis -- modified by parameterised volume
  //
  physiChamber = new G4PVParameterised(
                            "Chamber",       // their name
                            logicChamber,    // their logical volume
                            physiTracker,    // Mother physical volume
			    kZAxis,          // Are placed along this axis 
                            NbOfChambers,    // Number of chambers
                            chamberParam);   // The parametrisation

  G4cout << "There are " << NbOfChambers << " chambers in the tracker region. "
         << "The chambers are " << ChamberWidth/mm << " mm of " 
         << ChamberMater->GetName() << "\n The distance between chamber is "
	 << ChamberSpacing/cm << " cm" << G4endl;
	 
  //------------------------------------------------ 
  // Sensitive detectors
  //------------------------------------------------ 

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String trackerChamberSDname = "ExN02/TrackerChamberSD";
  ExN02TrackerSD* aTrackerSD = new ExN02TrackerSD( trackerChamberSDname );
  SDman->AddNewDetector( aTrackerSD );
  logicChamber->SetSensitiveDetector( aTrackerSD );

//--------- Visualization attributes -------------------------------

  G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  logicWorld  ->SetVisAttributes(BoxVisAtt);  
  logicTarget ->SetVisAttributes(BoxVisAtt);
  logicTracker->SetVisAttributes(BoxVisAtt);
  
  G4VisAttributes* ChamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  logicChamber->SetVisAttributes(ChamberVisAtt);
  
//--------- example of User Limits -------------------------------

  //below is an example of how to set tracking constraints in a given
  //logical volume(see also in N02PhysicsList how to setup the process
  //G4UserSpecialCuts).  
  //set a max Step length in the tracker region
  ///G4double maxStep = 0.5*ChamberWidth, maxLength = 2*fTrackerLength;
  ///G4double maxTime = 0.1*ns, minEkin = 10*MeV;
  ///logicTracker->SetUserLimits(new G4UserLimits(maxStep,maxLength,maxTime,
  ///                                              minEkin));
  
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
 
void ExN02DetectorConstruction::setTargetMaterial(G4String materialName)
{
  // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialName);  
  if (pttoMaterial)
     {TargetMater = pttoMaterial;
      logicTarget->SetMaterial(pttoMaterial); 
      G4cout << "\n----> The target is " << fTargetLength/cm << " cm of "
             << materialName << G4endl;
     }             
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ExN02DetectorConstruction::setChamberMaterial(G4String materialName)
{
  // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialName);  
  if (pttoMaterial)
     {ChamberMater = pttoMaterial;
      logicChamber->SetMaterial(pttoMaterial); 
      G4cout << "\n----> The chambers are " << ChamberWidth/cm << " cm of "
             << materialName << G4endl;
     }             
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..... 
 
void ExN02DetectorConstruction::SetMagField(G4double fieldValue)
{
  fpMagField->SetFieldValue(fieldValue);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..... 
