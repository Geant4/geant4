// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em2DetectorConstruction.cc,v 1.2 1999-12-15 14:49:00 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em2DetectorConstruction.hh"
#include "Em2DetectorMessenger.hh"

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "G4RunManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UnitsTable.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em2DetectorConstruction::Em2DetectorConstruction()
:solidEcal(NULL) ,logicEcal(NULL) ,physiEcal(NULL),
 solidSlice(NULL),logicSlice(NULL),physiSlice(NULL),
 solidRing(NULL) ,logicRing(NULL) ,physiRing(NULL), 
 myMaterial(NULL),magField(NULL)  ,
 nLtot(20),dLradl(1.),nRtot(20),dRradl(0.25),
 EcalLength(0.),EcalRadius(0.)
{
  detectorMessenger = new Em2DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em2DetectorConstruction::~Em2DetectorConstruction()
{ delete detectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* Em2DetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em2DetectorConstruction::DefineMaterials()
{
  G4String name, symbol; 
  G4double a, z, density;
  G4int ncomponents, natoms;
  G4double fractionmass;

  //
  // define few Elements
  //
  
    a = 1.01*g/mole;
    G4Element* H = new G4Element(name="Hydrogen", symbol="H", z=1., a);
      
    a = 14.01*g/mole;
    G4Element* N = new G4Element(name="Nitrogen", symbol="N", z=7., a);
    
    a = 16.00*g/mole;
    G4Element* O = new G4Element(name="Oxygen"  , symbol="O", z=8., a);
    
    a = 72.59*g/mole;
    G4Element* Ge = new G4Element(name="Germanium", symbol="Ge",z=32., a);
         
    a = 183.84*g/mole;
    G4Element* W = new G4Element(name="Tungsten"  , symbol="W", z=74., a);
    
    a = 207.19*g/mole;
    G4Element* Pb = new G4Element(name="Lead"     , symbol="Pb",z=82., a);
    
    a = 208.98*g/mole;
    G4Element* Bi = new G4Element(name="Bismuth"  , symbol="Bi",z=83., a);
      
  //  
  // define materials
  //
       
  //Air
    density = 1.29*mg/cm3;
    G4Material* Air = new G4Material(name="Air", density, ncomponents=2);
    Air->AddElement(N, fractionmass=0.7);
    Air->AddElement(O, fractionmass=0.3);
    
  //H2O
    density = 1.00*g/cm3;
    G4Material* H2O = new G4Material(name="Water", density, ncomponents=2);
    H2O->AddElement(H, natoms=2);
    H2O->AddElement(O, natoms=1);
    
  //liquid argon
    a = 39.95*g/mole;
    density = 1.390*g/cm3;
    G4Material* lAr = new G4Material(name="lAr", z=18., a, density);
        
  //Al
    a = 26.98*g/mole;
    density = 2.7*g/cm3;
    G4Material* Al = new G4Material(name="Al", z=13., a, density);
        
  //Fe
    a = 55.85*g/mole;
    density = 7.87*g/cm3;
    G4Material* Fe = new G4Material(name="Fe", z=26., a, density);
    
  //BGO
    density = 7.10*g/cm3;
    G4Material* BGO = new G4Material(name="BGO", density, ncomponents=3);
    BGO->AddElement(O , natoms=12);
    BGO->AddElement(Ge, natoms= 3);
    BGO->AddElement(Bi, natoms= 4); 
      
  //PbWO4
    density = 8.28*g/cm3;
    G4Material* PbWO = new G4Material(name="PbWO4", density, ncomponents=3);
    PbWO->AddElement(O , natoms=4);
    PbWO->AddElement(Pb, natoms=1);
    PbWO->AddElement(W , natoms=1);
     
  //Pb
    density = 11.35*g/cm3;
    G4Material* pb = new G4Material(name="Lead", density, ncomponents=1);
    pb->AddElement(Pb, fractionmass=1.0);
     
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  //choose material
  myMaterial = PbWO;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
G4VPhysicalVolume* Em2DetectorConstruction::ConstructVolumes()
{
  G4double Radl = myMaterial->GetRadlen();

  G4double dL = dLradl*Radl, dR = dRradl*Radl;  
  EcalLength = nLtot*dL; EcalRadius = nRtot*dR;
  
  //   
  // Ecal
  //
  solidEcal = new G4Tubs("Ecal",0.,EcalRadius,0.5*EcalLength,0.,360*deg);
  logicEcal = new G4LogicalVolume( solidEcal,myMaterial,"Ecal",0,0,0);
  physiEcal = new G4PVPlacement(0,G4ThreeVector(),
                                "Ecal",logicEcal,0,false,0);
				
  // Ring
  solidRing = new G4Tubs("Ring",0.,dR,0.5*EcalLength,0.,360*deg);
  logicRing = new G4LogicalVolume(solidRing,myMaterial,"Ring",0,0,0);
  if (nRtot >1)
     physiRing = new G4PVReplica("Ring",logicRing,physiEcal,
                               kRho,nRtot,dR);
  else
     physiRing = new G4PVPlacement(0,G4ThreeVector(),"Ring",logicRing,
                                    physiEcal,false,0);
				                   				
  // Slice
  solidSlice = new G4Tubs("Slice",0.,dR,0.5*dL,0.,360*deg);
  logicSlice = new G4LogicalVolume(solidSlice,myMaterial,"Slice",0,0,0);
  if (nLtot >1)
     physiSlice = new G4PVReplica("Slice",logicSlice,physiRing,
                                   kZAxis,nLtot,dL);
  else
     physiSlice = new G4PVPlacement(0,G4ThreeVector(),"Slice",logicSlice,
                                    physiRing,false,0);               


  cout << "Absorber is " << G4BestUnit(EcalLength,"Length") 
       << " of " << myMaterial->GetName() << G4endl; 


  G4VisAttributes* VisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  VisAtt->SetVisibility(true);
  logicEcal->SetVisAttributes(VisAtt);
  
  //
  //always return the physical World
  //
  return physiEcal;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em2DetectorConstruction::SetMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial)
     {myMaterial = pttoMaterial;
      logicEcal->SetMaterial(myMaterial); 
     }             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em2DetectorConstruction::SetLBining(G4ThreeVector Value)
{
  nLtot = (G4int)Value(0); dLradl = Value(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em2DetectorConstruction::SetRBining(G4ThreeVector Value)
{
  nRtot = (G4int)Value(0); dRradl = Value(1);
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em2DetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr 
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    
  if(magField) delete magField;		//delete the existing magn field
  
  if(fieldValue!=0.)			// create a new one if non nul
  { magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));        
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);
  } else {
    magField = NULL;
    fieldMgr->SetDetectorField(magField);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void Em2DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
