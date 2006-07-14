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
// $Id: A01DetectorConstruction.cc,v 1.1 2006-07-14 14:43:12 asaim Exp $
// --------------------------------------------------------------
//

#include "A01DetectorConstruction.hh"

///////////////////////////////////#include "G4FieldManager.hh"
///////////////////////////////////#include "G4TransportationManager.hh"
///////////////////////////////////#include "A01Field.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSNofStep.hh"
#include "G4PSTrackLength.hh"
#include "G4PSNofSecondary.hh"
#include "G4PSEnergyDeposit.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

A01DetectorConstruction::A01DetectorConstruction()
: air(0), water(0), worldVisAtt(0), phantomVisAtt(0) 
{;}

A01DetectorConstruction::~A01DetectorConstruction() 
{
    delete worldVisAtt;
    delete phantomVisAtt;
}

G4VPhysicalVolume* A01DetectorConstruction::Construct() 
{
    ConstructMaterials();

    G4Box * worldSolid = new G4Box("worldBox", 1.*m, 1.*m, 1.*m);
    G4LogicalVolume * worldLogical
	= new G4LogicalVolume(worldSolid, air, "worldLogical", 0, 0, 0);
    G4VPhysicalVolume* worldPhysical
	= new G4PVPlacement(0, G4ThreeVector(), worldLogical, "worldPhysical",
			    0, 0, 0);
    // water phantom
    G4Box * phantomSolid = new G4Box("phantomBox", 60.*cm, 60.*cm, 60*cm);
    G4LogicalVolume * phantomLogical
	= new G4LogicalVolume(phantomSolid, water, "phantomLogical", 0, 0 ,0);
    new G4PVPlacement(0, G4ThreeVector(), phantomLogical, "phantomPhysical",
		      worldLogical, 0, 0);
    G4Box * layerSolid = new G4Box("PhantomLayerBox", 60.*cm, 60.*cm, 3.*cm);
    G4LogicalVolume * layerLogical
        = new G4LogicalVolume(layerSolid, water, "PhantomLayerLogical", 0, 0, 0);
    new G4PVReplica("PhantomLayerPhysical",layerLogical,phantomLogical,kZAxis,20,6.*cm);

    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4String SDname;
    G4MultiFunctionalDetector* aSD = new G4MultiFunctionalDetector(SDname="MassWorld");
    SDman->AddNewDetector(aSD);
    layerLogical->SetSensitiveDetector(aSD);

    aSD->RegisterPrimitive(new G4PSNofStep("NofStep"));
    aSD->RegisterPrimitive(new G4PSTrackLength("TrackLength"));
    aSD->RegisterPrimitive(new G4PSNofSecondary("NofSecondary"));
    aSD->RegisterPrimitive(new G4PSEnergyDeposit("EnergyDeposit"));

    worldVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    worldLogical->SetVisAttributes(worldVisAtt);

    phantomVisAtt = new G4VisAttributes(G4Colour(0.6,0.8,1.0));
    phantomLogical->SetVisAttributes(phantomVisAtt);
    layerLogical->SetVisAttributes(phantomVisAtt);

    return worldPhysical;
}

void A01DetectorConstruction::ConstructMaterials() {

    G4double a;
    G4double z;
    G4double density;
    G4double weightRatio;
    G4String name;
    G4String symbol;
    G4int nElem;
    G4int nAtoms;

    // elements for mixtures and compounds
    a = 1.01*g/mole;
    G4Element* elH = new G4Element(name="Hydrogen", symbol="H", z=1., a);
    a = 14.01*g/mole;
    G4Element* elN = new G4Element(name="Nitrogen", symbol="N", z=7., a);
    a = 16.00*g/mole;
    G4Element* elO = new G4Element(name="Oxigen", symbol="O", z=8., a);

    // Air
    density = 1.29*mg/cm3;
    air = new G4Material(name="Air", density, nElem=2);
    air->AddElement(elN, weightRatio=.7);
    air->AddElement(elO, weightRatio=.3);

    // Water
    density = 1.0*g/cm3;
    water = new G4Material(name="H2O", density, nElem=2);
    water->AddElement(elH, nAtoms=2);
    water->AddElement(elO, nAtoms=1);

    G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}


