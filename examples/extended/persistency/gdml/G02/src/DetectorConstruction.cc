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
// $Id: DetectorConstruction.cc,v 1.8 2009-04-24 15:47:13 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Class DetectorConstruction implementation
//
// ----------------------------------------------------------------------------

#include "DetectorConstruction.hh"

// Geant4 includes
//
#include "globals.hh"
#include "G4GeometryManager.hh"
#include "G4VisAttributes.hh"

// Materials
//
#include "G4Material.hh"

// Geometry includes
//
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"

// Reflected solids
//
#include "G4ReflectedSolid.hh"
#include "G4DisplacedSolid.hh"
#include "G4ReflectionFactory.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4Transform3D.hh"

// Assembly volumes
//
#include "G4AssemblyVolume.hh"

// Volume parameterisations
//
#include "ChamberParameterisation.hh"

// Messenger
//
#include "DetectorMessenger.hh"

// GDML parser include
//
#include "G4GDMLParser.hh"

// ----------------------------------------------------------------------------
//
// Constructor
//
DetectorConstruction::DetectorConstruction()
  : Air(0), Aluminum(0), Pb(0), Xenon(0)
{
  expHall_x=5.*m;
  
  fReadFile ="test.gdml";
  fWriteFile="wtest.gdml";
  fStepFile ="mbb";
  writingChoice=1;
 
  detectorMessenger = new DetectorMessenger( this );
}

// ----------------------------------------------------------------------------
//
// Destructor
//
DetectorConstruction::~DetectorConstruction()
{
  if(detectorMessenger) delete detectorMessenger;
}

// ----------------------------------------------------------------------------
//
// Constructs geometries and materials
//
G4VPhysicalVolume* DetectorConstruction::Construct()
{ 
  // Writing or Reading of Geometry using G4GDML

  G4VPhysicalVolume* fWorldPhysVol;

  if(writingChoice==0)
  {
    // **** LOOK HERE*** FOR READING GDML FILES
    //
    
    // ACTIVATING OVERLAP CHECK when read volumes are placed.
    // Can take long time in case of complex geometries
    //
    // parser.SetOverlapCheck(true);

    parser.Read(fReadFile);

    // READING GDML FILES OPTION: 2nd Boolean argument "Validate".
    // Flag to "false" disables check with the Schema when reading GDML file.
    // See the GDML Documentation for more information.
    //
    // parser.Read(fReadFile,false);
     
    // Prints the material information
    //
    G4cout << *(G4Material::GetMaterialTable() ) << G4endl;
         
    // Giving World Physical Volume from GDML Parser
    //
    fWorldPhysVol = parser.GetWorldVolume();     
  }
  else if(writingChoice==1)
  {
    // **** LOOK HERE*** FOR WRITING GDML FILES
    // Detector Construction and WRITING to GDML
    //
    ListOfMaterials();
    fWorldPhysVol = ConstructDetector();

    // OPTION: TO ADD MODULE AT DEPTH LEVEL ...
    //
    // Can be a integer or a pointer to the top Physical Volume:
    //
    // G4int depth=1;
    // parser.AddModule(depth);
     
    // OPTION: SETTING ADDITION OF POINTER TO NAME TO FALSE
    //
    // By default, written names in GDML consist of the given name with
    // appended the pointer reference to it, in order to make it unique.
    // Naming policy can be changed by using the following method, or
    // calling Write with additional Boolean argument to "false".
    // NOTE: you have to be sure not to have duplication of names in your
    //       Geometry Setup.
    // 
    // parser.SetAddPointerToName(false);
    //
    // or
    //
    // parser.Write(fWriteFile, fWorldPhysVol, false);
    
    // Writing Geometry to GDML File
    //
    parser.Write(fWriteFile, fWorldPhysVol);
     
    // OPTION: SPECIFYING THE SCHEMA LOCATION
    //
    // When writing GDML file the default the Schema Location from the
    // GDML web site will be used:
    // "http://cern.ch/service-spi/app/releases/GDML/GDML_2_10_0/src/GDMLSchema/gdml.xsd"
    //
    // NOTE: GDML Schema is distributed in Geant4 in the directory:
    //    $G4INSTALL/source/persistency/gdml/schema
    //
    // You can change the Schema path by adding a parameter to the Write
    // command, as follows:
    //
    // parser.Write(fWriteFile, fWorldPhysVol, "your-path-to-schema/gdml.xsd");
  }
  else   // Demonstration how to Read STEP files using GDML 
  {
     // Some printout...
     //
     ListOfMaterials();

     // Arbitrary values that should enclose any reasonable geometry
     //
     const G4double expHall_y = expHall_x/50.;
     const G4double expHall_z = expHall_x/50.;

     // Create the hall
     //
     G4Box * experimentalHallBox
       = new G4Box("ExpHallBox",expHall_x/50.,expHall_y,expHall_z);
     G4LogicalVolume * experimentalHallLV
       = new G4LogicalVolume(experimentalHallBox, Air,"ExpHallLV");
     fWorldPhysVol
       = new G4PVPlacement(0, G4ThreeVector(0.0,0.0,0.0),
                           experimentalHallLV, "ExpHallPhys", 0, false, 0);

     // DetectorConstruction via reading STEP File
     //
     G4LogicalVolume* LogicalVolST
       = parser.ParseST(fStepFile,Air,Aluminum);

     // Placement inside of the hall
     //
     new G4PVPlacement(0, G4ThreeVector(10.0,0.0,0.0), LogicalVolST,
                       "StepPhys", experimentalHallLV, false, 0);
  }

  // Set Visualization attributes to world
  //
  G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  fWorldPhysVol->GetLogicalVolume()->SetVisAttributes(BoxVisAtt);  

  return fWorldPhysVol;
}

// ----------------------------------------------------------------------------
//
// Utility to build and list necessary materials
//
void DetectorConstruction::ListOfMaterials()
{
  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density,temperature,pressure;
  G4double fractionmass;
  G4String name, symbol;
  G4int ncomponents;

  // Elements needed for the materials

  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", z=7., a);

  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxygen", symbol="O", z=8., a);
          
  a = 26.98*g/mole;
  G4Element* elAl = new G4Element(name="Aluminum", symbol="Al", z=13., a);
        
  // Print the Element information
  //
  G4cout << *(G4Element::GetElementTable()) << G4endl;

  // Air
  //
  density = 1.29*mg/cm3;
  Air = new G4Material(name="Air", density, ncomponents=2);
  Air->AddElement(elN, fractionmass=0.7);
  Air->AddElement(elO, fractionmass=0.3);

  // Aluminum
  //
  density = 2.70*g/cm3;
  Aluminum = new G4Material(name="Aluminum", density, ncomponents=1);
  Aluminum->AddElement(elAl, fractionmass=1.0);

  // Lead
  //
  Pb = new G4Material("Lead", z=82., a= 207.19*g/mole, density= 11.35*g/cm3);

  // Xenon gas
  //
  Xenon = new G4Material("XenonGas", z=54., a=131.29*g/mole,
                         density= 5.458*mg/cm3, kStateGas,
                         temperature= 293.15*kelvin, pressure= 1*atmosphere);

  // Prints the material information
  //
  G4cout << *(G4Material::GetMaterialTable() ) << G4endl;
}

// ----------------------------------------------------------------------------
//
// Detector Construction
//
// Detector consist from DetectorBox, Conrol Room and 4 SubDetectors
// SubDetectors1 and 2 show how to use Reflection Factory and Assembly
// SubDetectors 3 and 4 show how to use Parameterisation
//
G4VPhysicalVolume* DetectorConstruction::ConstructDetector()
{
  // Arbitary values that should enclose any reasonable geometry
  //
  const G4double expHall_y = expHall_x;
  const G4double expHall_z = expHall_x;

  // Create the hall
  //
  G4Box * experimentalHallBox =
      new G4Box("ExpHallBox", expHall_x, expHall_y, expHall_z);
  G4LogicalVolume * experimentalHallLV =
      new G4LogicalVolume(experimentalHallBox, Air, "ExpHallLV");
  G4PVPlacement * experimentalHallPhys =
      new G4PVPlacement(0, G4ThreeVector(0.0,0.0,0.0), experimentalHallLV,
                        "ExpHallPhys", 0, false, 0);

  // DetectorConstruction

  const G4double det_x = expHall_x*0.8;
  const G4double det_y = expHall_x*0.7;
  const G4double det_z = det_y;

  // Create the detector box
  //
  G4Box * detectorBox =
      new G4Box("detectorBox", det_x, det_y, det_z);
  G4LogicalVolume * detectorLV =
      new G4LogicalVolume(detectorBox, Air, "detLV");
  // G4PVPlacement * detectorPhys = 
      new G4PVPlacement(0, G4ThreeVector(0.0,0.0,0.0), detectorLV,
                        "detPhys", experimentalHallLV, false, 0);

  // Create the Control room box
  //
  const G4double room_x = expHall_x/20.;
  const G4double room_y = room_x;
  const G4double room_z = room_x;

  G4Box * roomBox =
      new G4Box("roomBox", room_x, room_y, room_z);
  G4LogicalVolume * roomLV =
      new G4LogicalVolume(roomBox, Air, "roomLV");
  // G4PVPlacement * roomPhys =
      new G4PVPlacement(0, G4ThreeVector(expHall_x-room_x-10.,0.0,0.0), roomLV,
                        "roomPhys", experimentalHallLV, false, 0);
  
  // SubDetector1
  //
  const G4double bigL=expHall_x/5.+50.;
  G4LogicalVolume* subDetectorLV1 = ConstructSubDetector1();
  // G4PVPlacement * detPhys1 =
      new G4PVPlacement(0, G4ThreeVector(bigL,0.0,0.0), subDetectorLV1,
                        "PhysSubDetector1", detectorLV, false, 0); 
        
  //
  // LOOK HERE FOR REFLECTIONS
  //   

  // SubDetector2
  //     
  G4Translate3D translation(-bigL, 0., 0.);
  G4RotationMatrix* rotD3 = new G4RotationMatrix();
  G4Transform3D rotation = G4Rotate3D(*rotD3);
  G4ReflectX3D  reflection;
  G4Transform3D transform = translation*rotation*reflection;

  // Place the reflected part using G4ReflectionFactory
  //
  G4ReflectionFactory::Instance()->Place(transform, "reflSubDetector",
                                         subDetectorLV1, detectorLV, false, 0);

  // SubDetector3
  //
  G4LogicalVolume* subDetectorLV3 = ConstructSubDetector2();
  // G4PVPlacement * detPhys3 =
      new G4PVPlacement(0, G4ThreeVector(0.0,bigL,0.0), subDetectorLV3,
                        "PhysSubDetectorFirst3", detectorLV, false, 0); 
         
  // SubDetector4, placement of parameterised chambers
  //
  G4LogicalVolume* subDetectorLV4 = ConstructSubDetector2();
  G4LogicalVolume* subChamberLV   = ConstructParametrisationChamber();
  // G4PVPlacement * detChamb =
      new G4PVPlacement(0, G4ThreeVector(0,0.0,0.0), subChamberLV,
                        "AssemblyPhys", subDetectorLV4, false, 0);
  // G4PVPlacement * detPhys4 =
      new G4PVPlacement(0, G4ThreeVector(0.0,-bigL,0.0), subDetectorLV4,
                        "PhysSubDetectorSecond3", detectorLV, false, 0);

  return experimentalHallPhys;
}

// ----------------------------------------------------------------------------
//
// SubDetector1
//
G4LogicalVolume* DetectorConstruction::ConstructSubDetector1()
{    
  const G4double sub_x = expHall_x/5.;
  const G4double sub_y = sub_x;

  // Create the hall
  //
  G4Tubs * subTub =
    new G4Tubs("subTub", 0., sub_x, sub_y, -90.*deg, 180*deg);
  G4LogicalVolume * subTubLV =
    new G4LogicalVolume(subTub, Pb, "tubLV");
  G4LogicalVolume *AssemblyLV = ConstructAssembly();
  // G4PVPlacement * detAss =
    new G4PVPlacement(0, G4ThreeVector(sub_x/3,0.0,0.0), AssemblyLV,
                      "AssemblyPhys", subTubLV, false, 0);
  return subTubLV;  
}  

// ----------------------------------------------------------------------------
//
// SubDetector2
//
G4LogicalVolume* DetectorConstruction::ConstructSubDetector2()
{    
  const G4double sub_x = expHall_x/10.;
  const G4double sub_y = sub_x*2.;
  const G4double sub_z = sub_x;

  // Create the hall
  //
  G4Box * detHallBox =
    new G4Box("detHallBox", sub_x, sub_y, sub_z);
  G4LogicalVolume * detHallLV =
    new G4LogicalVolume(detHallBox, Aluminum, "detHallLV");

  return detHallLV;  
} 

// ----------------------------------------------------------------------------
//
// Assembly
//
G4LogicalVolume* DetectorConstruction::ConstructAssembly()
{    
  const G4double big_x = expHall_x/17;
  const G4double big_y = big_x;
  const G4double big_z = big_x;

  // Create the Box
  //
  G4Box * OuterBox =
    new G4Box("OuterBox", big_x, big_y, big_z);
  G4LogicalVolume * OuterBoxLV =
    new G4LogicalVolume(OuterBox, Air, "OuterBoxLV");
  // G4PVPlacement * OuterBoxPhys = 
    new G4PVPlacement(0, G4ThreeVector(0.0,0.0,0.0), OuterBoxLV,
                      "OuterBoxPhys", 0, false, 0);

  // The aluminum object's logical volume
  //
  const G4double bigL=big_x/2.5;
  const G4double medL=big_x/8;
  const G4double smalL=big_x/12;

  G4Box * BigBox =
    new G4Box("BBox", bigL, bigL, bigL);
  G4LogicalVolume * BigBoxLV =
    new G4LogicalVolume(BigBox, Aluminum, "AlBigBoxLV");
  G4Box * MedBox =
    new G4Box("MBox", medL, medL, medL);
  G4LogicalVolume * MedBoxLV1 =
    new G4LogicalVolume(MedBox, Aluminum, "AlMedBoxLV1");
  G4Box * SmallBox =
    new G4Box("SBox", smalL, smalL, smalL);
  G4LogicalVolume * SmallBoxLV =
    new G4LogicalVolume(SmallBox, Aluminum, "AlSmaBoxLV");

  const G4double bigPlace=bigL+10.;
  const G4double medPlace=medL+10.;
  // G4PVPlacement * BigBoxPhys = 
    new G4PVPlacement(0, G4ThreeVector(bigPlace,0.0,0.0), BigBoxLV,
                      "AlPhysBig", OuterBoxLV, false, 0); 

  // Construction of Tub 
  //
  G4Tubs * BigTube =
    new G4Tubs("BTube",0,smalL,smalL,-pi/2.,pi);

  // Construction of Reflection of Tub
  //
  G4ReflectX3D  Xreflection;
  G4Translate3D translation(-bigPlace, 0., 0.);
  G4Transform3D transform =Xreflection;

  G4ReflectedSolid * ReflBig =
    new G4ReflectedSolid("Refll_Big", BigTube, transform);
  G4LogicalVolume * ReflBigLV =
    new G4LogicalVolume(ReflBig, Xenon, "ReflBigAl");
   new G4PVPlacement(0, G4ThreeVector(0.,0.0,0.0), ReflBigLV,
                      "AlPhysBigTube", SmallBoxLV, false, 0); 
  //
  // LOOK HERE FOR ASSEMBLY 
  //

  // create Assembly of Boxes and Tubs
  //
  G4AssemblyVolume* assembly = new G4AssemblyVolume();
  G4RotationMatrix* rot = new G4RotationMatrix();
  G4ThreeVector posBig(-bigPlace, 0, 0);
  G4ThreeVector posBig0(bigPlace/4, 0, 0);
  G4ThreeVector posMed(-medPlace, 0, 0);
  G4ThreeVector posMed0(medPlace, 0, 0);
  G4ThreeVector position(0., 0., 0.);

  // Add to Assembly the MediumBox1
  //
  assembly->AddPlacedVolume(MedBoxLV1, posMed0, rot);

  // Add to Assembly the Small Box
  //
  assembly->AddPlacedVolume(SmallBoxLV, posMed, rot);
  
  // Place the Assembly
  //
  assembly->MakeImprint(BigBoxLV, posBig0, rot, 0);

  //
  // LOOK HERE FOR ASSEMBLY with REFLECTION 
  //   

  G4Translate3D translation1(-bigPlace, 0., 0.);
  G4RotationMatrix* rotD3 = new G4RotationMatrix();
  G4Transform3D rotation = G4Rotate3D(*rotD3);
  G4ReflectX3D  reflection;
  G4Transform3D transform1 = translation1*rotation*reflection;

  assembly->MakeImprint(OuterBoxLV, transform1, 0, 0);

  return OuterBoxLV;  
}

// ----------------------------------------------------------------------------
//
// Parameterised Chamber
//
G4LogicalVolume* DetectorConstruction::ConstructParametrisationChamber()
{    
  const G4double chamber_x = expHall_x/12.;
  const G4double chamber_y = chamber_x;
  const G4double chamber_z = chamber_x;

  // Create the hall
  //
  G4Box * paramChamberBox =
    new G4Box("ChamberBox", chamber_x, chamber_y, chamber_z);
  G4LogicalVolume * paramChamberLV =
    new G4LogicalVolume(paramChamberBox, Air, "ChamberLV");

  // Parametrisation Chamber (taken from N02 novice example)
  //
  G4int NbOfChambers = 5;
  G4double ChamberWidth = 2*cm;
  G4double ChamberSpacing = 8*cm;
  G4double fTrackerLength = (NbOfChambers+1)*ChamberSpacing; // Full length
  G4double trackerSize = 0.5*fTrackerLength;

  // An example of parameterised volume
  // dummy values for G4Box -- modified by parameterised volume
  //
  G4Box *solidChamber =
    new G4Box("chamber", 10*cm, 10*cm, 1*cm); 
  G4LogicalVolume* logicChamber =
    new G4LogicalVolume(solidChamber, Aluminum, "Chamber", 0, 0, 0);
    
  G4double firstPosition = -trackerSize + 0.5*ChamberWidth;
  G4double firstLength = fTrackerLength/10;
  G4double lastLength  = fTrackerLength;

  G4VPVParameterisation* chamberParam =
    new ChamberParameterisation( NbOfChambers,          // NoChambers 
                                 firstPosition,         // Z of center of first 
                                 ChamberSpacing,        // Z spacing of centers
                                 ChamberWidth,          // Width Chamber 
                                 firstLength,           // lengthInitial 
                                 lastLength);           // lengthFinal
  // G4VPhysicalVolume* physiChamber = 
    new G4PVParameterised( "Chamber",          // their name
                           logicChamber,       // their logical volume
                           paramChamberLV,     // mother logical volume      
                           kZAxis,             // Are placed along this axis 
                           NbOfChambers,       // Number of chambers
                           chamberParam);      // The parametrisation
  return paramChamberLV;
} 

// ----------------------------------------------------------------------------
//
// SetReadFile
//
void DetectorConstruction::SetReadFile( const G4String& File )
{
  fReadFile=File;
  writingChoice=0;
}

// ----------------------------------------------------------------------------
//
// SetWriteFile
//
void DetectorConstruction::SetWriteFile( const G4String& File )
{
  fWriteFile=File;
  writingChoice=1;
}

// ----------------------------------------------------------------------------
//
// SetStepFile
//
void DetectorConstruction::SetStepFile( const G4String& File )
{
  fStepFile=File;
  writingChoice=3;
}
