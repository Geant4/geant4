// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TstVADetectorConstruction.cc,v 1.3 2001-02-07 17:31:01 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------

#include "TstVADetectorConstruction.hh"
#include "TstVADetectorMessenger.hh"

#include "g4std/strstream"

#include "G4ios.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4TransportationManager.hh"
#include "G4GeometryManager.hh"
#include "G4StateManager.hh"
#include "G4UImanager.hh"
#include "G4TransportationManager.hh"
#include "G4AssemblyVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"

TstVADetectorConstruction::TstVADetectorConstruction()
:worldVol(0),Air(0),Al(0),Pb(0),selectedMaterial(0),detectorChoice(0),plateLV(0)
{
  classicDetector.caloLV      = 0;
  classicDetector.PVs.clear()    ;
  assemblyDetector            = 0;
  ConstructClassic();
  materialChoice              = "Air";
  detectorMessenger           = new TstVADetectorMessenger(this);
}

TstVADetectorConstruction::~TstVADetectorConstruction()
{
  // R.I.P. messenger
  delete detectorMessenger;
  CleanClassic();
  CleanAssembly(); 
}

G4VPhysicalVolume* TstVADetectorConstruction::Construct()
{
  if( worldVol == 0 )
  {
    switch(detectorChoice)
    { 
      case 1:
        ConstructAssembly(); 
        break;
      default:
        ConstructClassic();
    }
  }
  return worldVol;
}

void TstVADetectorConstruction::SwitchDetector()
{
  CleanClassic();
  CleanAssembly();
  switch(detectorChoice)
  { 
    case 1:
      ConstructAssembly(); 
      break;
    default:
      ConstructClassic(); 
  }

  // Notify run manager that the new geometry has been built
  G4RunManager::GetRunManager()->DefineWorldVolume( worldVol );
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  // Let the navigator to know about the new top of the new geometry
  G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->SetWorldVolume(worldVol);
  // We have to tweak the navigator's state. By the following dummy call we ensure that navigator's
  // state is reset properly
  G4ThreeVector center(0,0,0);
  G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->LocateGlobalPointAndSetup(center,0,false);
}

void TstVADetectorConstruction::SelectDetector(G4String val)
{
  if(val=="assembly") 
  { detectorChoice = 1; }
  else
  { detectorChoice = 0; }
  G4cout << "Now Detector is " << val << G4endl;
}

void TstVADetectorConstruction::SelectMaterial(G4String val)
{
  materialChoice = val;
  SelectMaterialPointer();
  G4cout << "World volume is now made of " << materialChoice << G4endl;
}

void TstVADetectorConstruction::SelectMaterialPointer()
{
//--------- Material definition ---------

  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;

  if(!Air)
  {
    a = 14.01*g/mole;
    G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a);
    a = 16.00*g/mole;
    G4Element* elO = new G4Element(name="Oxigen", symbol="O", iz=8., a);
    density = 1.29e-03*g/cm3;
    Air = new G4Material(name="Air", density, nel=2);
    Air->AddElement(elN, .7);
    Air->AddElement(elO, .3);
  }

  if(!Al)
  {
    a = 26.98*g/mole;
    density = 2.7*g/cm3;
    Al = new G4Material(name="Aluminium", z=13., a, density);
  }

  if(!Pb)
  {
    a = 207.19*g/mole;
    density = 11.35*g/cm3;
    Pb = new G4Material(name="Lead", z=82., a, density);
  }

  if(materialChoice=="Air")
  { selectedMaterial = Air; }
  else if(materialChoice=="Al")
  { selectedMaterial = Al; }
  else
  { selectedMaterial = Pb; }

  G4LogicalVolume* worldLV = worldVol->GetLogicalVolume();
  if(worldLV)
  { worldLV->SetMaterial(selectedMaterial); }
}

const double worldX              = 2000*mm;
const double worldY              = 2000*mm;
const double worldZ              = 2000*mm;

const double caloX               = 1600*mm;
const double caloY               = 1600*mm;
const double caloZ               =  200*mm;

const double plateX              =  700*mm;
const double plateY              =  700*mm;
const double plateZ              =  100*mm;

const unsigned int layers        =    5;

const double firstCaloPos        =  500*mm;
const double caloCaloOffset      =   50*mm;
const double plateCaloOffset     =    1*mm;
const double platePlateOffset    =    2*mm;

void TstVADetectorConstruction::ConstructClassic()
{
  if( worldVol == 0 )
  {
    // Define world volume
    G4Box*           WorldBox  = new G4Box( "WBox", worldX/2., worldY/2., worldZ/2. );
    G4LogicalVolume* worldLV   = new G4LogicalVolume( WorldBox, selectedMaterial, "WLog", 0, 0, 0);
    worldVol                   = new G4PVPlacement(0, G4ThreeVector(), "WPhys", worldLV, 0, false, 0);

    // We need to this here to overcome the chicken-egg problem of proper initialization of the
    // world volume material
    if( selectedMaterial == 0 )
    {
      SelectMaterialPointer();
    }
    
    // Define a calorimeter layer
    G4Box*           CaloBox   = new G4Box( "CaloBox", caloX/2., caloY/2., caloZ/2. );
    classicDetector.caloLV     = new G4LogicalVolume( CaloBox, Air, "CaloLV", 0, 0, 0 );

    // Define a calorimeter plate
    G4Box*           PlateBox  = new G4Box( "PlateBox", plateX/2., plateY/2., plateZ/2. );
    plateLV                    = new G4LogicalVolume( PlateBox, Pb, "PlateLV", 0, 0, 0 );

    // Fill layer with plates
    G4VPhysicalVolume* platePV = new G4PVPlacement(  (G4RotationMatrix*)0
                                                    ,G4ThreeVector( caloX/4. 
                                                                   ,caloY/4. 
                                                                   ,0.)
                                                    ,plateLV
                                                    ,"plate_pv_0"
                                                    ,classicDetector.caloLV
                                                    ,false, 0 );

    // Remember the instance so we can clean up properly later
    classicDetector.PVs.push_back( platePV );

    platePV                    = new G4PVPlacement(  (G4RotationMatrix*)0
                                                    ,G4ThreeVector( -1*caloX/4. 
                                                                   ,caloY/4.
                                                                   ,0.)
                                                    ,plateLV
                                                    ,"plate_pv_1"
                                                    ,classicDetector.caloLV
                                                    ,false, 1 );

    // Remember the instance so we can clean up properly later
    classicDetector.PVs.push_back( platePV );

    platePV                    = new G4PVPlacement(  (G4RotationMatrix*)0
                                                    ,G4ThreeVector( -1*caloX/4. 
                                                                   ,-1*caloY/4.
                                                                   ,0.)
                                                    ,plateLV
                                                    ,"plate_pv_2"
                                                    ,classicDetector.caloLV
                                                    ,false, 2 );

    // Remember the instance so we can clean up properly later
    classicDetector.PVs.push_back( platePV );

    platePV                    = new G4PVPlacement(  (G4RotationMatrix*)0
                                                    ,G4ThreeVector( caloX/4.
                                                                   ,-1*caloY/4. 
                                                                   ,0.)
                                                    ,plateLV
                                                    ,"plate_pv_3"
                                                    ,classicDetector.caloLV
                                                    ,false, 3 );

    // Remember the instance so we can clean up properly later
    classicDetector.PVs.push_back( platePV );

    // Create layers of quazi calorimeter
    for( unsigned int i = 0; i < layers; i++ )
    {
      G4std::strstream pvName;

      pvName << "CaloPV_" << i << G4std::ends;
      
      // Place each layer
      G4VPhysicalVolume* caloPV = new G4PVPlacement( 0
                                                    ,G4ThreeVector
                                                      ( 0
                                                       ,0
                                                       ,i*( caloZ + caloCaloOffset) - firstCaloPos
                                                      )
                                                    ,pvName.str(), classicDetector.caloLV, worldVol,
                                                     false, i );

      // Remember the instance so we can clean up properly later
      classicDetector.PVs.push_back( caloPV );
    }

#ifdef G4DEBUG    
    G4cout << "PVs created: " << classicDetector.PVs.size() << G4endl;
#endif
  }
}

void TstVADetectorConstruction::ConstructAssembly()
{
  if( worldVol == 0 )
  {
    // Define world volume
    G4Box*           WorldBox  = new G4Box( "WBox", worldX/2., worldY/2., worldZ/2. );
    G4LogicalVolume* worldLV   = new G4LogicalVolume( WorldBox, selectedMaterial, "WLog", 0, 0, 0);
    worldVol                   = new G4PVPlacement(0, G4ThreeVector(), "WPhys", worldLV, 0, false, 0);

    // We need to this here to overcome the chicken-egg problem of proper initialization of the
    // world volume material
    if( selectedMaterial == 0 )
    {
      SelectMaterialPointer();
    }

    // Define a calorimeter plate
    G4Box*           PlateBox  = new G4Box( "PlateBox", plateX/2., plateY/2., plateZ/2. );
    plateLV                    = new G4LogicalVolume( PlateBox, Pb, "PlateLV", 0, 0, 0 );
    
    // Define one calorimeter layer as one assembly volume
    assemblyDetector           = new G4AssemblyVolume();

    // Rotation and translation of a plate inside the assembly
    G4RotationMatrix        Ra;
    G4ThreeVector           Ta;

    // Rotation of the assembly inside the world
    G4RotationMatrix        Rm;

    // Fill the assembly by the plates  
    Ta.setX( caloX/4. );  Ta.setY( caloY/4. );  Ta.setZ( 0. );
    assemblyDetector->AddPlacedVolume( plateLV, Ta, &Ra );
    
    Ta.setX( -1*caloX/4. );  Ta.setY( caloY/4. );  Ta.setZ( 0. );
    assemblyDetector->AddPlacedVolume( plateLV, Ta, &Ra );
    
    Ta.setX( -1*caloX/4. );  Ta.setY( -1*caloY/4. );  Ta.setZ( 0. );
    assemblyDetector->AddPlacedVolume( plateLV, Ta, &Ra );

    Ta.setX( caloX/4. );  Ta.setY( -1*caloY/4. );  Ta.setZ( 0. );
    assemblyDetector->AddPlacedVolume( plateLV, Ta, &Ra );
    
    // Now instantiate the layers of calorimeter
    for( unsigned int i = 0; i < layers; i++ )
    {
      // Translation of the assembly inside the world
      G4ThreeVector Tm( 0,0,i*(caloZ + caloCaloOffset) - firstCaloPos );
      assemblyDetector->MakeImprint( worldLV, Tm, &Rm );
    }
  }
}

void TstVADetectorConstruction::CleanClassic()
{
  // First free the memory occupied by physical volumes
  for( unsigned int i = 0; i < classicDetector.PVs.size(); i++ )
  {
    G4VPhysicalVolume* toDie  = classicDetector.PVs[i];
    if( toDie != 0 )
    {
      // Clean up the rotation matrix if any
      G4RotationMatrix* rmToDie = toDie->GetRotation();
      if( rmToDie != 0 )
      {
        delete rmToDie;
      }
      delete toDie;
    }
  }
  
  classicDetector.PVs.clear();
  
  // Now free the memory of logical volume objects
  if( classicDetector.caloLV != 0 )
  {
    G4VSolid* solToDie = classicDetector.caloLV->GetSolid();
    if( solToDie != 0 )
    {
      delete solToDie;
    }
    delete classicDetector.caloLV;
    classicDetector.caloLV = 0;
  }

  if( plateLV != 0 )
  {
    G4VSolid* solToDie = plateLV->GetSolid();
    if( solToDie != 0 )
    {
      delete solToDie;
    }
    delete plateLV;
    plateLV = 0;
  }

  // Finally R.I.P. world instance
  if( worldVol != 0 )
  {
    G4LogicalVolume* worldLV = worldVol->GetLogicalVolume();
    if( worldLV != 0 )
    {
       G4VSolid* solToDie = worldLV->GetSolid();
       if( solToDie != 0 )
       {
         delete solToDie;
       }
       delete worldLV;
    }
  
    delete worldVol;
    worldVol = 0;
  }
}

void TstVADetectorConstruction::CleanAssembly()
{
  // Clean-up of assembly is simple :-)
  if( assemblyDetector != 0 )
  {
    delete assemblyDetector;
    assemblyDetector = 0;
  }

  // Clean the plates logical volume
  if( plateLV != 0 )
  {
    G4VSolid* solToDie = plateLV->GetSolid();
    if( solToDie != 0 )
    {
      delete solToDie;
    }
    delete plateLV;
    plateLV = 0;
  }
  
  // Finally R.I.P. world instance
  if( worldVol != 0 )
  {
    G4LogicalVolume* worldLV = worldVol->GetLogicalVolume();
    if( worldLV != 0 )
    {
       G4VSolid* solToDie = worldLV->GetSolid();
       if( solToDie != 0 )
       {
         delete solToDie;
       }
       delete worldLV;
    }
  
    delete worldVol;
    worldVol = 0;
  }
}

