// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3toG4DetectorConstruction.cc,v 1.1 2000-07-24 11:23:49 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//--------------------------------------------------------------------------
// G3toG4DetectorConstruction. Most the work is Done in
// G4BuildGeom, which returns a G4LogicalVolume*, a pointer to the
// top-level logiical volume in the detector defined by the call List file
// inFile
//--------------------------------------------------------------------------

#include "G4ios.hh"
#include "G3toG4DetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Material.hh"
#include "G4Box.hh"

G3toG4DetectorConstruction::G3toG4DetectorConstruction(G4String inFile){
    _inFile = inFile;
    G4cout << "Instantiated G3toG4DetectorConstruction using call list file \""
           << _inFile << "\"" << G4endl;
}

G3toG4DetectorConstruction::~G3toG4DetectorConstruction(){
  //    G4cout << "Deleted G3toG4DetectorConstruction..." << G4endl;
}

G4VPhysicalVolume*
G3toG4DetectorConstruction::Construct(){
  _lv = G4BuildGeom(_inFile);
  //_lv = SimpleConstruct();
  if (_lv != 0) {
    _pv = new G4PVPlacement(0, G4ThreeVector(), _lv, _lv->GetName(), 0,
			    false, 0);
    G4cout << "Top-level G3toG4 logical volume " << _lv->GetName() << " "
	   << *(_lv -> GetVisAttributes()) << G4endl;
  } else 
    G4cerr << "creation of logical mother failed !!!" << G4endl;
  return _pv;
}
G4LogicalVolume*
G3toG4DetectorConstruction::SimpleConstruct(){
  G4String name, symbol;             //a=mass of a mole;
  G4double a, z, density, fractionmass; //z=mean number of protons;  
  G4int ncomponents, iz, n;          //iz=number of protons  in an isotope; 
  // n=number of nucleons in an isotope;

  a = 14.01*g/mole;
  G4Element* N  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

  a = 16.00*g/mole;
  G4Element* O  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

  //
  // define a material from elements.   case 2: mixture by fractional mass
  //

  density = 1.290*mg/cm3;
  G4Material* Air = new G4Material(name="Air"  , density, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);
  G4VSolid* Mother = new G4Box("TestMother",		//its name
				100*cm, 100*cm, 100*cm); //its size

  G4VSolid* Daughter = new G4Box("TestDaughter", 50*cm, 20*cm, 10*cm);
    			     
  G4LogicalVolume* logicMother = new G4LogicalVolume(Mother,	//its solid
						     Air,	//its material
						     "LTestMother");//its name
    				       
  G4LogicalVolume* logicDaughter = new G4LogicalVolume(Daughter, //its solid
						       Air,	//its material
						       "LTestDaughter"); 

  G4VPhysicalVolume* physiDaughter = new G4PVPlacement(0,
						       G4ThreeVector(),
						       logicDaughter,
						       "PTestDaughter",
						       logicMother, 
						       false,0);
  //                                        
  // Visualization attributes
  //

  logicMother->SetVisAttributes (G4VisAttributes::Invisible);
  G4VisAttributes* DaughterVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  DaughterVisAtt->SetVisibility(true);
  logicDaughter->SetVisAttributes(DaughterVisAtt);
  return logicMother;
}













