// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T08DetectorConstruction.cc,v 1.1 1999-01-08 16:35:19 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "T08DetectorConstruction.hh"
#include "T08DetectorMessenger.hh"
#include "T08ChamberParameterisation.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"
#include "T08MagneticField.hh"
#include "G4PVParameterised.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
///#include "G4UserLimits.hh"
#include "G4SDManager.hh"
#include "T08TrackerSD.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

T08DetectorConstruction::T08DetectorConstruction()
:solidWorld(NULL),logicWorld(NULL),physiWorld(NULL),
 solidTracker(NULL),logicTracker(NULL),physiTracker(NULL), 
 solidTarget(NULL),logicTarget(NULL),physiTarget(NULL), 
 solidChamber(NULL),logicChamber(NULL),physiChamber(NULL), 
 myMaterial(NULL), fDetectorLength(0.),fWorldLength(0.)
{
  detectorMessenger = new T08DetectorMessenger(this);
  fpMagField = new T08MagneticField(); // G4FieldManager();
  // magFieldManager = new G4FieldManager();
}

T08DetectorConstruction::~T08DetectorConstruction()
{ 
  delete detectorMessenger;
  delete fpMagField;           
  // delete magFieldManager;
}

G4VPhysicalVolume* T08DetectorConstruction::Construct()
{
//--------- Material definition ---------

  G4double a, iz, z, density;
  G4String name, symbol;
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
    
  //Vacuum
    a = 1.*g/mole;
    density = 1.e-20*g/cm3;
    G4Material* Vacuum = new G4Material(name="Vacuum",z=1.,a,density);  

//--------- Sizes of the principal geometrical components (solids)  ---------
  
  fDetectorLength = 1250.*cm ;       // Full length of the detector
  fTrackerLength = 500. * cm;       // Full length of the Tracker
  fTargetLength = 5.0 * cm;         // Full length of the Target
  
  G4double detectorSize = 0.5*fDetectorLength;
  G4double trackerSize = fTrackerLength * 0.5;   // Half length of the Tracker
  G4double targetSize  = fTargetLength  * 0.5;   // Half length of the Target
  
  
//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
  myMaterial = Aluminium;
  
  //------------------------------ 
  // World
  //------------------------------ 
  fWorldLength= 1.2 * fDetectorLength; 
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
  solidTarget = new G4Box("TargetSolid",targetSize,targetSize,targetSize);
  logicTarget = new G4LogicalVolume(solidTarget  ,myMaterial,"TargetLV",0,0,0);
  physiTarget = new G4PVPlacement(0,               // no rotation
				  G4ThreeVector(), // at (0,0,0)
				  "TargetPV",      // its name
				  logicTarget,     // its logical volume
				  physiWorld,      // its mother  volume
				  false,           // no boolean operations
				  0);              // no particular field 

  G4cout << "Target is " << 2*targetSize/cm << " cm of Al " << endl;

  //------------------------------ 
  // Tracker
  //------------------------------ 
  G4ThreeVector positionTracker;
  solidTracker = new G4Box("Tracker",trackerSize,trackerSize,trackerSize);
  logicTracker = new G4LogicalVolume(solidTracker , Air, "Tracker",0,0,0);
  positionTracker= G4ThreeVector(0,0,trackerSize+targetSize); 
  
  physiTracker = new G4PVPlacement(0,              // no rotation
				  positionTracker, // at (0,0,0)
				  "Tracker",       // its name
				  logicTracker,    // its logical volume
				  physiWorld,      // its mother  volume
				  false,           // no boolean operations
				  0);              // no particular field 

  // Ensure that the tracker fits inside the detector volume
  assert( 2.*trackerSize+targetSize <= detectorSize );

  //------------------------------ 
  // Tracker segments
  //------------------------------ 
  // An example of Parameterised volumes
  // dummy values for G4Box -- modified by parameterised volume
  solidChamber
    = new G4Box("chamberBox", 100*cm, 100*cm, 10*cm); 
                      // chamberWidth, chamberWidth, chamberDepth);
                 
  G4LogicalVolume * trackerChamberLV
    = new G4LogicalVolume(solidChamber, Aluminium, "trackerChamberLV",0,0,0);
  G4VPVParameterisation * chamberParam
    = new T08ChamberParameterisation(  
			      6,            //  NoChambers, 
			   -240.*cm,        //  Z of center of first 
			     80.*cm,        //  Z spacing of centers
			     20.*cm,        //  Width Chamber, 
			     50.*cm,        //  lengthInitial, 
			     2.*trackerSize); //  lengthFinal
  // dummy value : kXAxis -- modified by parameterised volume
  G4VPhysicalVolume *trackerChamber_phys
    = new G4PVParameterised("TrackerChamber_parameterisedPV",
                            trackerChamberLV, // Its logical volume
                            physiTracker,     // Mother physical volume
			    kZAxis,           // Are placed along this axis 
                            6,                // Number of chambers
                            chamberParam);    // The parametrisation
  G4VisAttributes* trackerChamber_logVisAtt
    = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
  trackerChamberLV->SetVisAttributes(trackerChamber_logVisAtt);

  //------------------------------------------------ 
  // Sensitive detectors
  //------------------------------------------------ 

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String trackerChamberSDname = "T08/TrackerChamberSD";
  T08TrackerSD* aTrackerSD = new T08TrackerSD( trackerChamberSDname );
  SDman->AddNewDetector( aTrackerSD );
  trackerChamberLV->SetSensitiveDetector( aTrackerSD );

//--------- Visualization attributes -------------------------------

  logicWorld->SetVisAttributes (G4VisAttributes::Invisible);

  G4VisAttributes * simpleBoxVisAtt
    = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  logicTracker->SetVisAttributes(simpleBoxVisAtt);

  //set a max Step length in the absorber
  ///G4double maxStep = fDetectorLength/0.7, maxLength=fDetectorLength/0.7;
  ///G4double maxTime = 0.1*ns, minEkin = 10*MeV;
  ///logicAbsorber->SetUserLimits(new G4UserLimits(maxStep,maxLength,maxTime,
  ///                                              minEkin));
  
  return physiWorld;
}

#if 0
void T08DetectorConstruction::setMaterial(G4String materialChoice)
{
   const G4MaterialTable* 
                       theMaterialTable = G4Material::GetMaterialTable();
   
   G4Material* pttoMaterial;
   for (G4int J=0 ; J<theMaterialTable->length() ; J++)
    { pttoMaterial = (*theMaterialTable)(J);     
      if(pttoMaterial->GetName() == materialChoice)
         {myMaterial = pttoMaterial;
          logicTracker->SetMaterial(pttoMaterial); 
          G4cout << "Absorber is " << fDetectorLength/cm << " cm of " << materialChoice 
               << endl;}             
    }

}
#endif 

void T08DetectorConstruction::SetDetectorLength(G4double length)
{
  G4GeometryManager::GetInstance()->OpenGeometry();
  
  fDetectorLength = length; 
  fWorldLength = 1.2*fDetectorLength;
  
  G4double halfSize = 0.5 * fDetectorLength; 
  G4double halfWorldLength = 0.5*fWorldLength;
  
  solidWorld -> SetXHalfLength(halfWorldLength);
  solidWorld -> SetYHalfLength(halfWorldLength);  
  solidWorld -> SetZHalfLength(halfWorldLength); 
  
  solidTracker -> SetXHalfLength(halfSize);
  solidTracker -> SetYHalfLength(halfSize);  
  solidTracker -> SetZHalfLength(halfSize);
  
  G4cout << "Absorber is " << fDetectorLength/cm << " cm of " << myMaterial->GetName()
       << endl;     
}  

void T08DetectorConstruction::SetMagField(G4double fieldValue)
{
  fpMagField->SetFieldValue(fieldValue);
}
