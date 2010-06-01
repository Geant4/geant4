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

#include "Tst01DetectorConstruction.hh"

#include "G4ios.hh"
#include "Tst01DetectorMessenger.hh"

#include <sstream>

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4Torus.hh"

#include "G4BooleanSolid.hh"
#include "G4DisplacedSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4Transform3D.hh"
#include "G4VoxelLimits.hh"

#include "G4LogicalVolume.hh"
#include "G4AssemblyVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4TransportationManager.hh"
#include "G4GeometryManager.hh"
#include "G4StateManager.hh"
#include "G4UImanager.hh"
#include "G4RunManager.hh"
#include "G4TransportationManager.hh"

//////////////////////////////////////////////////////////////////
//
// Constructor/Destructor

Tst01DetectorConstruction::Tst01DetectorConstruction()
 : simpleBoxLog(0), simpleBoxDetector(0),
   honeycombDetector(0), fWorldPhysVol(0),
   fTestCSG(0), fTestLog(0), fTestVol(0),
   fDisPb(0), fPb1(0), fPb2(0), fPb3(0), fSphere(0),
   fTestBool(0), fTestBoolLog(0), fTestBoolVol(0),
   fTestD1Log(0), fTestD1Vol(0), Air(0), Al(0), Pb(0),
   selectedMaterial(0), detectorChoice(0), fChoiceCSG(0), fChoiceBool(0),
   AssemblyDetectorLog(0), AssemblyDetector(0), AssemblyCalo(0),
   AssemblyCellLog(0), AssemblyDetector2(0), AssemblyDetector3(0)
{
  wSize = 2000.*cm;

  detectorMessenger = new Tst01DetectorMessenger(this);
  materialChoice = "Pb";
}

Tst01DetectorConstruction::~Tst01DetectorConstruction()
{
  delete detectorMessenger;

  // Clean up Assembly setup
  //
  if( AssemblyCalo != 0 )  { delete AssemblyCalo; }
}

//////////////////////////////////////////////////////////////////////
//
// 

G4VPhysicalVolume* Tst01DetectorConstruction::Construct()
{
  if(!fWorldPhysVol)
  {
    ConstructDetectors();
  }

  switch(detectorChoice)
  { 
    case 1:
    {
      fWorldPhysVol = honeycombDetector; 
      break;
    }
    case 2:
    {
      fWorldPhysVol = AssemblyDetector; 
      break;
    }
    case 3:
    {
      fWorldPhysVol = AssemblyDetector2; 
      break;
    }
    case 4:
    {
      fWorldPhysVol = AssemblyDetector3; 
      break;
    }
    default:
    {
      fWorldPhysVol = simpleBoxDetector;
    }
  }
  return fWorldPhysVol;
}

////////////////////////////////////////////////////////////////////////
//
//

void Tst01DetectorConstruction::SwitchDetector()
{
  if(!fWorldPhysVol)
  {
    ConstructDetectors();
  }

  switch(detectorChoice)
  { 
    case 1:
    {
      fWorldPhysVol = honeycombDetector ;
      break;
    }
    case 2:
    {
      fWorldPhysVol = AssemblyDetector ;
      break;
    }
    case 3:
    {
      fWorldPhysVol = AssemblyDetector2 ;
      break;
    }
    case 4:
    {
      fWorldPhysVol = AssemblyDetector3 ;
      break;
    }
    default:
    {
      fWorldPhysVol = simpleBoxDetector ;
    }
  }
  G4RunManager::GetRunManager()->DefineWorldVolume(fWorldPhysVol);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

/////////////////////////////////////////////////////////////////////
//
//

void Tst01DetectorConstruction::SelectDetector(const G4String& val)
{
  if(val == "Honeycomb")
  {
    detectorChoice = 1 ;
  }
  else if(val == "Assembly") 
  {
    detectorChoice = 2 ;
  }
  else if(val == "Assembly2") 
  {
    detectorChoice = 3 ;
  }
  else if(val == "Assembly3") 
  {
    detectorChoice = 4 ;
  }
  else
  {
    detectorChoice = 0 ;
  }

  G4cout << G4endl;
  G4cout << "           -------------------------" << G4endl ;
  G4cout << "           Now Detector is " << val << G4endl ;
  G4cout << "           -------------------------" << G4endl ;
  G4cout << G4endl;
}
/////////////////////////////////////////////////////////////////////
//
//

void Tst01DetectorConstruction::SelectCSG(const G4String& name)
{
  if     ( name == "Tubs"   ) { fChoiceCSG = 1 ; }
  else if( name == "Cons"   ) { fChoiceCSG = 2 ; }
  else if( name == "Sphere" ) { fChoiceCSG = 3 ; }
  else                        { fChoiceCSG = 0 ; } // default or error in name
 
  G4cout << "---> CGS solid is " << name << G4endl;
}

////////////////////////////////////////////////////////////////////////
//
//

void Tst01DetectorConstruction::SwitchCSG()
{
  SelectMaterialPointer();

  if( fTestCSG ) { delete fTestCSG; }
  
  switch(fChoiceCSG)
  { 
    case 1 :
    {
      fTestCSG = new G4Tubs("testCSG", 20*cm, 50*cm, 50*cm, 0, 2*pi) ;
      break ;
    }
    case 2 :
    {
      fTestCSG = new G4Cons("testCSG",10*cm,30*cm, 20*cm,50*cm, 50*cm, 0, 2*pi) ;
      break ;
    }
    case 3 :
    {
      fTestCSG = new G4Sphere("testCSG",20*cm,50*cm, 0,2*pi, 0,2*pi) ;
      break ;
    }
    default:
    {
      fTestCSG = new G4Box("testCSG", 20*cm, 50*cm, 50*cm ) ;
    }
  }
  if( fTestLog ) { delete fTestLog; }
  fTestLog = new G4LogicalVolume(fTestCSG,selectedMaterial,"testLog",0,0,0) ;

  if( fTestVol ) { delete fTestVol; }
  fTestVol = new G4PVPlacement(0, G4ThreeVector(), "testVol", fTestLog, 
                               fWorldPhysVol, false, 0);
  
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

/////////////////////////////////////////////////////////////////////
//
//

void Tst01DetectorConstruction::SelectBoolean(const G4String& name)
{
  if(      name == "Union"   )     { fChoiceBool = 1 ; }
  else if( name == "Subtraction" ) { fChoiceBool = 2 ; }
  else                             { fChoiceBool = 0 ; }
 
  G4cout << "---> Boolean operation is " << name << G4endl;
}

////////////////////////////////////////////////////////////////////////
//
//

void Tst01DetectorConstruction::SwitchBoolean()
{
  SelectMaterialPointer();
  if( !fWorldPhysVol ) { ConstructDetectors(); }

  if( !fDisPb )
  {
    G4RotationMatrix identity ;

    fPb1 = new G4Box("b1",50*cm,50*cm,50*cm) ;
    fPb2 = new G4Box("b2",10*cm,10*cm,60*cm) ;
    fPb3 = new G4Box("b4",50*cm,50*cm, 5*cm) ;

    fDisPb = new G4DisplacedSolid("disPb", fPb3, &identity,
                                  G4ThreeVector(0,0,60*cm)) ;

    fSphere = new G4Sphere("sphere", 0, 10*cm, 0 , 2*pi, 0, pi);
  }

  G4ThreeVector putDaughter;

  if( fTestBool ) { delete fTestBool; }

  switch(fChoiceBool)
  { 
    case 1 :
    {
      fTestBool = new G4UnionSolid("b1UniondisPb", fPb1, fDisPb) ;
      break ;
    }
    case 2 :
    {
      fTestBool = new G4SubtractionSolid("b1SubtractB2", fPb1, fPb2) ;
      putDaughter = G4ThreeVector(0.,30*cm,0.);
      break ;
    }
    default:
    {
      fTestBool = new G4IntersectionSolid("b1IntersectB2", fPb1, fPb2) ;
    }
  }

  if( fTestBoolLog ) { delete fTestBoolLog; }
  fTestBoolLog = new G4LogicalVolume(fTestBool, selectedMaterial,
                                     "testBoolLog", 0, 0, 0) ;

  if( fTestBoolVol ) { delete fTestBoolVol; }
  fTestBoolVol = new G4PVPlacement(0, G4ThreeVector(), "testBoolVol",
                                   fTestBoolLog, fWorldPhysVol, false, 0);

  // daughters
  
  if( fTestD1Log ) { delete fTestD1Log; }
  fTestD1Log = new G4LogicalVolume(fSphere, selectedMaterial,
                                   "testD1Log", 0, 0, 0) ;

  if( fTestD1Vol ) { delete fTestD1Vol; }
  fTestD1Vol = new G4PVPlacement(0, putDaughter, "testD1Vol",
                                 fTestD1Log, fTestBoolVol, false, 0);

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//////////////////////////////////////////////////////////////////////
//
//

void Tst01DetectorConstruction::SelectMaterial(const G4String& val)
{
  materialChoice = val;
  SelectMaterialPointer();
  G4cout << G4endl;
  G4cout << "---> Daughter CSG/Boolean will be made of "
         << materialChoice << G4endl;
  G4cout << G4endl;
}

///////////////////////////////////////////////////////////////////////
//
//

void Tst01DetectorConstruction::SelectMaterialPointer()
{

//
// Material definition 
//

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

  if(simpleBoxLog)
  { simpleBoxLog->SetMaterial(selectedMaterial); }
}

/////////////////////////////////////////////////////////////////////////
//
//

void Tst01DetectorConstruction::ConstructDetectors()
{

//
//   simpleBoxDetector 
//

  G4Box * mySimpleBox = new G4Box("SBox",wSize, wSize, wSize);
  simpleBoxLog = new G4LogicalVolume( mySimpleBox,Air,
                                      "SLog",0,0,0);
  simpleBoxDetector = new G4PVPlacement(0,G4ThreeVector(),"SPhys",
                                        simpleBoxLog,0,false,0);

//
//  honeycombDetector
//

  size_t i,j;

  G4Box *myWorldBox= new G4Box("WBox",wSize, wSize, wSize);
  G4Box *myCalBox = new G4Box("CBox",1500*cm, 1500*cm, 1000*cm);
  G4Tubs *myTargetTube 
     = new G4Tubs("TTube",0*cm, 22.5*cm, 1000*cm, 0.*deg, 360.*deg);

  G4LogicalVolume *myWorldLog = new G4LogicalVolume(myWorldBox,Air,
                                                    "WLog", 0, 0, 0);
  G4LogicalVolume *myCalLog   = new G4LogicalVolume(myCalBox,Al,
                                                    "CLog", 0, 0, 0);
  G4LogicalVolume *myTargetLog= new G4LogicalVolume(myTargetTube,Pb,
                                                    "TLog", 0, 0, 0);

  honeycombDetector = new G4PVPlacement(0,G4ThreeVector(),"WPhys",
                                        myWorldLog,0,false,0);

  G4PVPlacement *myCalPhys=new G4PVPlacement(0,G4ThreeVector(),
                                               "CalPhys",
                                               myCalLog,
                                               honeycombDetector,
                                               false,0);
  G4double offset=22.5*cm, xTlate, yTlate;
  G4int copyNo=0;
  for (j=1;j<=25;j++)
  {
    yTlate = -1000.0*cm - 40.0*cm + j*80.0*cm;
    for (i=1;i<=50;i++)
    {
      std::stringstream tName1;
      tName1 << "TPhysA_" << j << "_" << i  << std::ends;
      xTlate = -1000.0*cm - 20.0*cm + i*45.0*cm - offset;
      // G4PVPlacement *myTargetPhysA =
         new G4PVPlacement(0,G4ThreeVector(xTlate,yTlate,0*cm),
                           tName1.str().c_str(),myTargetLog,myCalPhys,
                           false,copyNo++);
    }
  }
  for (j=1;j<=26;j++)
  {
    yTlate = -1000.0*cm - 80.0*cm + j*80.0*cm;
    for (i=1;i<=50;i++)
    {
      std::stringstream tName2;
      tName2 << "TPhysB_" << j << "_" << i  << std::ends;
      xTlate = -1000.0*cm - 20.0*cm + i*45.0*cm;
      // G4PVPlacement *myTargetPhysB =
         new G4PVPlacement(0,G4ThreeVector(xTlate,yTlate,0*cm),
                           tName2.str().c_str(),myTargetLog,myCalPhys,
                           false,copyNo++);
    }
  }

//
//  AssemblyDetector
//

  const double worldX              = wSize;
  const double worldY              = wSize;
  const double worldZ              = wSize;

  const double caloX               = 1600*mm;
  const double caloY               = 1600*mm;
  const double caloZ               =  200*mm;

  const double plateX              =  700*mm;
  const double plateY              =  700*mm;
  const double plateZ              =  100*mm;

  const unsigned int layers        =    5;

  const double firstCaloPos        =  500*mm;
  const double caloCaloOffset      =   50*mm;

  // Define world volume for Assembly detector
  //
  G4Box* AssemblyWorldBox =
    new G4Box( "AssemblyWorldBox", worldX, worldY, worldZ );
  G4LogicalVolume* AssemblyWorldLog =
    new G4LogicalVolume( AssemblyWorldBox, selectedMaterial,
                        "AssemblyWorldLog", 0, 0, 0);
  // G4VPhysicalVolume* AssemblyWorld    =
    new G4PVPlacement(0, G4ThreeVector(), "AssemblyWorld",
                      AssemblyWorldLog, 0, false, 0);

  // The Assembly detector
  //
  G4Box* AssemblyBox   =
    new G4Box( "AssemblyBox", worldX/2, worldY/2, worldZ/2 );
  AssemblyDetectorLog  =
    new G4LogicalVolume( AssemblyBox, selectedMaterial,
                        "AssemblyDetectorLog", 0, 0, 0);
  AssemblyDetector     =
    new G4PVPlacement(0, G4ThreeVector(), AssemblyDetectorLog,
                      "AssemblyDetector",
                      AssemblyWorldLog, false, 0);

  // Define a calorimeter plate
  G4Box* AssemblyCellBox =
    new G4Box( "AssemblyCellBox", plateX/2., plateY/2., plateZ/2. );
  AssemblyCellLog        =
    new G4LogicalVolume( AssemblyCellBox, Pb, "AssemblyCellLog", 0, 0, 0 );
  
  // Define one calorimeter layer as one assembly volume
  AssemblyCalo = new G4AssemblyVolume();

  // Rotation and translation of a plate inside the assembly
  G4RotationMatrix        Ra;
  G4ThreeVector           Ta;

  // Rotation of the assembly inside the world
  G4RotationMatrix        Rm;

  // Fill the assembly by the plates  
  Ta.setX( caloX/4. );  Ta.setY( caloY/4. );  Ta.setZ( 0. );
  AssemblyCalo->AddPlacedVolume( AssemblyCellLog, Ta, &Ra );
  
  Ta.setX( -1*caloX/4. );  Ta.setY( caloY/4. );  Ta.setZ( 0. );
  AssemblyCalo->AddPlacedVolume( AssemblyCellLog, Ta, &Ra );
  
  Ta.setX( -1*caloX/4. );  Ta.setY( -1*caloY/4. );  Ta.setZ( 0. );
  AssemblyCalo->AddPlacedVolume( AssemblyCellLog, Ta, &Ra );

  Ta.setX( caloX/4. );  Ta.setY( -1*caloY/4. );  Ta.setZ( 0. );
  AssemblyCalo->AddPlacedVolume( AssemblyCellLog, Ta, &Ra );
  
  // Now instantiate the layers of calorimeter
  //
  for( i = 0; i < layers; i++ )
  {
    // Translation of the assembly inside the world
    G4ThreeVector Tm( 0,0,i*(caloZ + caloCaloOffset) - firstCaloPos );
    AssemblyCalo->MakeImprint( AssemblyDetectorLog, Tm, &Rm, 0 );
  }
  G4cout << "Total imprinted physical volumes for AssemblyDetector-1: "
         << AssemblyCalo->TotalImprintedVolumes() << G4endl;

  //
  //  AssemblyDetector2 (assembly of assemblies)
  //
  
  G4VSolid* worldS
    = new G4Box("worldS", wSize, wSize, wSize);
  G4LogicalVolume* worldV2
    = new G4LogicalVolume(worldS, Air, "worldV2");
  AssemblyDetector2
    = new G4PVPlacement(0, G4ThreeVector(), worldV2,
                        "world2", 0, false, 0); 
   
  // Make the elementary assembly of the whole structure
  // Use mother volume instead of assembly as Geant4 does not 
  // support assebly of assemblies

  G4int ntooth = 5;
  G4double xplate = 25.*cm;
  G4double yplate = 50.*cm;
  G4double xtooth = 10.*cm;
  G4double ytooth = 0.5*yplate/ntooth;
  G4double dshift = 2.*xplate + xtooth;
  G4double xt,yt;

  G4AssemblyVolume* tplate = new G4AssemblyVolume();

  // plate volume
  G4Box* plateS
    = new G4Box("plateS", xplate, yplate, 1.*cm);
  G4LogicalVolume* plateV
    = new G4LogicalVolume(plateS, Al, "PLATE");
   
  // tooth volume
  G4Box* toothS
    = new G4Box("toothS", xtooth, ytooth, 1.*cm);
  G4LogicalVolume* toothV
    = new G4LogicalVolume(toothS, Al, "TOOTH");
   
  // compose assembly 
  G4ThreeVector pos0(0.,0., 0.);
  tplate->AddPlacedVolume(plateV, pos0, 0);
  for (G4int i=0; i<ntooth; i++) {
    xt = xplate+xtooth;
    yt = -yplate + (4*i+1)*ytooth;
    G4ThreeVector pos1(xt, yt, 0);
    tplate->AddPlacedVolume(toothV, pos1, 0);

    xt = -xplate-xtooth;
    yt = -yplate + (4*i+3)*ytooth;
    G4ThreeVector pos2(xt, yt, 0);
    tplate->AddPlacedVolume(toothV, pos2, 0);
  }  
  
  G4RotationMatrix* rot1 = new G4RotationMatrix();
  rot1->rotateX(90.*deg);
  G4RotationMatrix *rot;
  G4AssemblyVolume* cell = new G4AssemblyVolume();

  // Make a hexagone cell out of 6 toothplates. These can zip togeather
  // without generating overlaps (they are self-contained)
  for (G4int i2=0; i2<6; i2++) {
    G4double phi =  60.*i2 * deg;
    G4double xp = dshift*std::sin(phi);
    G4double yp = -dshift*std::cos(phi);
    rot = new G4RotationMatrix(*rot1);
    rot->rotateZ(phi); 
    G4ThreeVector pos(xp, yp,0.);
    cell->AddPlacedAssembly(tplate, pos, rot);
  }   

  // Make a row as an assembly of cells, then combine rows in a honeycomb
  // structure. This again works without any need to define rows as
  // "overlapping"
  G4AssemblyVolume* row = new G4AssemblyVolume();
  G4int ncells = 5;
  for (G4int i3=0; i3<ncells; i3++) {
    G4double ycell = (2*i3+1)*(dshift+10.*cm);
    G4ThreeVector pos1(0., ycell, 0.);
    row->AddPlacedAssembly(cell, pos1, 0);
    G4ThreeVector pos2(0., -ycell, 0.);
    row->AddPlacedAssembly(cell, pos2, 0);
  }

  G4double dxrow = 3.*(dshift+10.*cm)*std::tan(30.*deg);
  G4double dyrow = dshift+10.*cm;
  G4int nrows = 5;
  for (G4int i4=0; i4<nrows; i4++) {
    G4double xrow = 0.5*(2*i4+1)*dxrow;
    G4double yrow = 0.5*dyrow;
    if ((i4%2)==0) yrow = -yrow;
    G4ThreeVector pos1(xrow, yrow, 0.);
    row->MakeImprint(worldV2, pos1, 0, 0);
    G4ThreeVector pos2(-xrow, -yrow, 0.);
    row->MakeImprint(worldV2, pos2, 0, 0);
  }
  G4cout << "Total imprinted physical volumes for AssemblyDetector-2: "
         << row->TotalImprintedVolumes() << G4endl;

  //
  //  AssemblyDetector3 (assembly with reflections)
  //
  
  // World
  //
  G4LogicalVolume* worldV3
    = new G4LogicalVolume(worldS, Air, "worldV3");
  AssemblyDetector3
    = new G4PVPlacement(0, G4ThreeVector(), worldV3, "world3", 0, false, 0); 
   
  // Assembly 
  //
  G4AssemblyVolume* assembly = new G4AssemblyVolume();


  // Volumes
  //
  G4VSolid* consS
    = new G4Cons("consS", 10.*cm, 40.*cm, 20.*cm, 60.*cm, 50*cm, 0., 360.*deg);
  G4LogicalVolume* consV
    = new G4LogicalVolume(consS, Al, "CONS");

  // Place volume in assembly
  //
  HepGeom::Transform3D transform1
   = HepGeom::Translate3D( 110.*cm,0., 0.)
   * HepGeom::RotateY3D( 90.*deg);
  assembly->AddPlacedVolume(consV, transform1);

  HepGeom::Transform3D transform2
   = HepGeom::ReflectX3D()
   * HepGeom::Translate3D( 110.*cm,0., 0.)
   * HepGeom::RotateY3D( 90.*deg);
  assembly->AddPlacedVolume(consV, transform2);

  HepGeom::Transform3D transform3
   = HepGeom::Translate3D( 0., 110.*cm, 0.)
   * HepGeom::RotateX3D(-90.*deg);
  assembly->AddPlacedVolume(consV, transform3);

  HepGeom::Transform3D transform4
   = HepGeom::ReflectY3D()
   * HepGeom::Translate3D( 0., 110.*cm, 0.)
   * HepGeom::RotateX3D(-90.*deg);
  assembly->AddPlacedVolume(consV, transform4);

  // Make imprint
  //
  G4RotationMatrix* rotv = 0;
  G4ThreeVector posv;
  assembly->MakeImprint(worldV3, posv, rotv, 0);

  G4cout << "Total imprinted physical volumes for AssemblyDetector-3: "
         << assembly->TotalImprintedVolumes() << G4endl;

  std::vector<G4VPhysicalVolume*>::iterator viter = assembly->GetVolumesIterator();
  G4cout << "Imprinted volumes: ";
  for (size_t i=0; i<assembly->TotalImprintedVolumes(); i++)
  {
     G4cout << (*viter)->GetName() << " - "; viter++;
  }
  G4cout << G4endl;

//
// Visualization attributes 
// 

  G4VisAttributes * simpleBoxVisAtt
    = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  simpleBoxLog->SetVisAttributes(simpleBoxVisAtt);

  G4VisAttributes * experimentalHallVisAtt
    = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  experimentalHallVisAtt->SetVisibility(false);
  myWorldLog->SetVisAttributes(experimentalHallVisAtt);

  G4VisAttributes * calorimeterBoxVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  calorimeterBoxVisAtt->SetForceWireframe(true);
  calorimeterBoxVisAtt->SetVisibility(true);
  myCalLog->SetVisAttributes(calorimeterBoxVisAtt);

  G4VisAttributes * calorimeterTubeVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.5,0.5));
  calorimeterTubeVisAtt->SetForceWireframe(true);
  calorimeterTubeVisAtt->SetVisibility(true);
  myTargetLog->SetVisAttributes(calorimeterTubeVisAtt);

}

