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
// $Id: Tst10DetectorConstruction.cc,v 1.16 2010-06-03 19:56:34 tnikitin Exp $
// ------------------------------------------------------------
//  GEANT 4 class header file 
//
//      This is a version for maximum particle set
//  History
//        first version              09  Sept. 1998 by S.Magni
// ------------------------------------------------------------

#include "Tst10DetectorConstruction.hh"
#include "Tst10DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4RotationMatrix.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Hype.hh"
#include "G4Para.hh"
#include "G4Paraboloid.hh"
#include "G4Torus.hh"
#include "G4Trd.hh"
#include "G4GenericTrap.hh"
#include "G4Polyhedra.hh"
#include "G4Polycone.hh"
#include "G4TwistedBox.hh"
#include "G4TwistedTrap.hh"
#include "G4TwistedTrd.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4TwoVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4TransportationManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4GeometryManager.hh"
#include "G4StateManager.hh"
#include "G4UImanager.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"

Tst10DetectorConstruction::Tst10DetectorConstruction()
  : aVolume(0), PhysicalVolume(0), Water(0),
    Water1(0), aSurface(0), fHallSize(1*m)
{
  detectorMessenger = new Tst10DetectorMessenger (this);
}

Tst10DetectorConstruction::~Tst10DetectorConstruction()
{
}

G4bool Tst10DetectorConstruction::CleanGeometry()
{
  if (PhysicalVolume)
  {
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::Clean();
    G4LogicalVolumeStore::Clean();
    G4SolidStore::Clean();

    return true;
  }
  else
  {
    return false;
  }
}

void Tst10DetectorConstruction::SwitchDetector()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(PhysicalVolume);
}

G4VPhysicalVolume*
Tst10DetectorConstruction::SelectDetector( const G4String& val )
{
  //------------------- A Volume ----------------------

  if (val == "Sphere")
    aVolume = new G4Sphere ( "aSphere", 8.0*cm, 10.0*cm, 
                             0.0*deg, 360.0*deg,0.0*deg, 130.0*deg);
  else if (val == "Orb")
   aVolume = new G4Orb ( "aOrb", 10.0*cm );
  else if (val == "Box")          
    aVolume = new G4Box ( "aBox", 10*cm, 10*cm, 10*cm );
  else if (val == "Cone")        
    aVolume = new G4Cons ( "aCone", 2*cm, 6*cm, 8*cm, 14*cm,
                           10*cm, 10*deg, 300*deg ); 
  else if (val == "Tube")
    aVolume = new G4Tubs ( "aTube", 5*cm, 10*cm, 7*cm, 70*deg, 100*deg);
  else if (val == "Hype")
    aVolume = new G4Hype ("aHype", 10*cm, 20*cm, 0*deg, 360*deg, 10*cm );
  else if (val == "Torus")
    aVolume = new G4Torus ("aTorus", 10*cm, 15*cm, 20*cm, 0*deg, 60*deg);
  else if (val == "Para")
    aVolume = new G4Para ("aPara", 8*cm, 10*cm, 12*cm, 30*deg, 45*deg, 60*deg);
  else if (val == "Paraboloid")
    aVolume = new G4Paraboloid ("aParaboloid", 8*cm, 1*cm, 12*cm);
  else if (val == "Trd")
    aVolume = new G4Trd ("aTrd", 8*cm, 10*cm, 7*cm, 9*cm, 10*cm);
  else if (val == "GenericTrap" ){
   std::vector<G4TwoVector> vertices;
   vertices.push_back(G4TwoVector( -4.5*cm, -4.5*cm));
   vertices.push_back(G4TwoVector( -4.5*cm,  4.5*cm));
   vertices.push_back(G4TwoVector(  4.5*cm,  4.5*cm));
   vertices.push_back(G4TwoVector(  4.5*cm, -4.5*cm));
   vertices.push_back(G4TwoVector( -3.5*cm, -3.5*cm));
   vertices.push_back(G4TwoVector( -3.5*cm,  3.5*cm));
   vertices.push_back(G4TwoVector(  3.5*cm,  3.5*cm));
   vertices.push_back(G4TwoVector(  3.5*cm, -2.5*cm));     
   aVolume = new G4GenericTrap("aGenTrd",14.*cm,vertices);
  }
  else if (val == "Polyhedra")
  {
    G4double zPlane[2] = {-10.*cm, 10.*cm };
    G4double rInner[2] = { 0.*cm,  0.*cm };
    G4double rOuter[2] = { 10.*cm, 10.*cm };
    aVolume = new G4Polyhedra("aPhedra",0*deg,360*deg,6,2,zPlane,rInner,rOuter);
  }
  else if (val == "Polycone")
  {
    G4double zPlane[2] = {-10.*cm, 10.*cm };
    G4double rInner[2] = { 2.*cm,  6.*cm };
    G4double rOuter[2] = { 10.*cm, 12.*cm };
    aVolume = new G4Polycone("aPcone",0*deg,360*deg,2,zPlane,rInner,rOuter);
  }
  else if (val == "TwistedBox")
  {
    aVolume = new G4TwistedBox("aTwistedBox",40*deg,5*cm,10*cm,15*cm);
  }
  else if (val == "TwistedTrd")
  {
    aVolume = new G4TwistedTrd("aTwistedTrd",5*cm,10*cm,8*cm,15*cm,18*cm,20*deg);
  }
  else if (val == "TwistedTrap")
  {
    aVolume = new G4TwistedTrap("aTwistedTrap",40*deg,5*cm,10*cm,8*cm,15*cm);
  }
  else if ( val == "TwistedTrap2") 
  {
    aVolume = new G4TwistedTrap("aTwistedTrap2",
				   20*deg,    // twist angle
				   80*cm,         // half z length
				   10*deg,      // direction between end planes
				   40*deg,        // defined by polar and azimutal angles.
				   8*cm,        // half y length at -pDz
				   11*cm,        // half x length at -pDz,-pDy
				   16*cm,        // half x length at -pDz,+pDy
				   8*cm,        // half y length at +pDz
				   11*cm,         // half x length at +pDz,-pDy
				   16*cm,        // half x length at +pDz,+pDy
				   -50*deg        // tilt angle at +pDz
				   ) ;
  }
  else
  {
    G4Exception("Tst10DetectorConstruction::SelectDetector() - Invalid shape!");
  }

  G4Box * Hall
          = new G4Box("Hall", fHallSize, fHallSize, fHallSize );
  G4LogicalVolume * Hall_log
          = new G4LogicalVolume (Hall, Water, "Hall_L", 0,0,0);
  PhysicalVolume
          = new G4PVPlacement(0,G4ThreeVector(),"Hall_P",Hall_log,0,false,0);

   
  G4LogicalVolume * aVolume_log 
    = new G4LogicalVolume(aVolume, Water1, "aVolume_L", 0,0,0);

   G4RotationMatrix *mat=new G4RotationMatrix();
  mat->rotateX(30.*deg); 

  G4VPhysicalVolume * aVolume_phys1
    = new G4PVPlacement(0,G4ThreeVector(50*cm, 0*cm, 0*cm),val, 
                        aVolume_log, PhysicalVolume, false, 0);
  G4VPhysicalVolume * aVolume_phys2
    = new G4PVPlacement(mat,G4ThreeVector(-50*cm, 0*cm, 0*cm),val, 
                        aVolume_log, PhysicalVolume, false, 0);
  G4VPhysicalVolume * aVolume_phys3
    = new G4PVPlacement(mat,G4ThreeVector(0*cm, 50*cm, 0*cm),val, 
                        aVolume_log, PhysicalVolume, false, 0);
  G4VPhysicalVolume * aVolume_phys4
    = new G4PVPlacement(mat,G4ThreeVector(0*cm, -50*cm, 0*cm),val, 
                        aVolume_log, PhysicalVolume, false, 0);
  G4VPhysicalVolume * aVolume_phys5
    = new G4PVPlacement(0,G4ThreeVector(0*cm, 0*cm, 50*cm),val, 
                        aVolume_log, PhysicalVolume, false, 0);
  G4VPhysicalVolume * aVolume_phys6
    = new G4PVPlacement(0,G4ThreeVector(0*cm, 0*cm, -50*cm),val, 
                    aVolume_log, PhysicalVolume, false, 0);

  // ------------ Surfaces definition ------------------

  G4LogicalBorderSurface* BorderSurfaces[12];   
  BorderSurfaces[0] = new G4LogicalBorderSurface("VolumeSurface",
                               PhysicalVolume,
                               aVolume_phys1,
                               aSurface);
  BorderSurfaces[1] = new G4LogicalBorderSurface("VolumeSurface",
                               PhysicalVolume,
                               aVolume_phys2,
                               aSurface);
  BorderSurfaces[2] = new G4LogicalBorderSurface("VolumeSurface",
                               PhysicalVolume,
                               aVolume_phys3,
                               aSurface);
  BorderSurfaces[3] = new G4LogicalBorderSurface("VolumeSurface",
                               PhysicalVolume,
                               aVolume_phys4,
                               aSurface);
  BorderSurfaces[4] = new G4LogicalBorderSurface("VolumeSurface",
                               PhysicalVolume,
                               aVolume_phys5,
                               aSurface);
  BorderSurfaces[5] = new G4LogicalBorderSurface("VolumeSurface",
                               PhysicalVolume,
                               aVolume_phys6,
                               aSurface);
  BorderSurfaces[6] = new G4LogicalBorderSurface("VolumeSurface",
                               aVolume_phys1,
                               PhysicalVolume,
                               aSurface);
  BorderSurfaces[7] = new G4LogicalBorderSurface("VolumeSurface",
                               aVolume_phys2,
                               PhysicalVolume,
                               aSurface);
  BorderSurfaces[8] = new G4LogicalBorderSurface("VolumeSurface",
                               aVolume_phys3,
                               PhysicalVolume,
                               aSurface);
  BorderSurfaces[9] = new G4LogicalBorderSurface("VolumeSurface",
                               aVolume_phys4,
                               PhysicalVolume,
                               aSurface);
  BorderSurfaces[10] = new G4LogicalBorderSurface("VolumeSurface",
                               aVolume_phys5,
                               PhysicalVolume,
                               aSurface);
  BorderSurfaces[11] = new G4LogicalBorderSurface("VolumeSurface",
                               aVolume_phys6,
                               PhysicalVolume,
                               aSurface);

  G4cout << "\n You selected a " << val << " detector" << G4endl;

  return PhysicalVolume;
}

void Tst10DetectorConstruction::SetMaterial()
{
  G4String name, symbol;
  G4double density = 1.00*g/cm3;
  G4double a, iz;
  Water = new G4Material(name="Water", density, 2);
  Water1 = new G4Material(name="Water1", density, 2);
  a = 1*g/mole;
  G4Element* elH = new G4Element(name="Hydrogen", symbol="H", iz=1., a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxygen", symbol="O", iz=8., a);

  Water->AddElement(elH, .66);
  Water->AddElement(elO, .34);
  Water1->AddElement(elH, .66);
  Water1->AddElement(elO, .34);

  const G4int NUMENTRIES = 5;
  G4double RINDEX_WATER [NUMENTRIES];
  G4double RINDEX_WATER1 [NUMENTRIES];
  G4double REFLECTIVITY [NUMENTRIES];
  G4double EFFICIENCY [NUMENTRIES];
  
  for (int i=0; i<NUMENTRIES; i++) {
    RINDEX_WATER1[i]=5.0;
    RINDEX_WATER[i]=1.33;
    REFLECTIVITY[i]=0.9;
    EFFICIENCY[i]=1.0;
  }  
  G4double PHENERGY[NUMENTRIES] =
            { 0.01, 1.0, 2.0, 3.0, 4.0};
  G4MaterialPropertiesTable *WaterMPT = new G4MaterialPropertiesTable();
  G4MaterialPropertiesTable *WaterMPT1 = new G4MaterialPropertiesTable();
  WaterMPT->AddProperty("RINDEX", PHENERGY, RINDEX_WATER, NUMENTRIES);
  WaterMPT1->AddProperty("RINDEX", PHENERGY, RINDEX_WATER1, NUMENTRIES);
  Water->SetMaterialPropertiesTable(WaterMPT);
  Water1->SetMaterialPropertiesTable(WaterMPT1);

  aSurface = new G4OpticalSurface ( "aSurface" );
  aSurface->SetType(dielectric_metal);
  aSurface->SetFinish(polishedfrontpainted);
  aSurface->SetModel(glisur);  
  G4MaterialPropertiesTable* SurfaceMPT = new G4MaterialPropertiesTable();
  SurfaceMPT->AddProperty("REFLECTIVITY", PHENERGY, REFLECTIVITY, NUMENTRIES);
  SurfaceMPT->AddProperty("EFFICIENCY", PHENERGY, EFFICIENCY, NUMENTRIES);
  aSurface->SetMaterialPropertiesTable ( SurfaceMPT );
}

G4VPhysicalVolume* Tst10DetectorConstruction::Construct()
{
  SetMaterial();

  //-------------------Hall ----------------------------------
  
  return SelectDetector ("Sphere");
}
