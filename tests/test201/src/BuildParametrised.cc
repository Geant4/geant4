// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: BuildParametrised.cc,v 1.1 1999-01-08 16:35:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Based on test code from Hans-Peter Wellisch

#include "BuildParametrised.hh"

#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"

#include "G4PVParameterised.hh"
#include "ParametrisedBox.hh"

G4VPhysicalVolume* BuildParametrised ()
{
  G4RotationMatrix* theNull = new G4RotationMatrix ();
  G4ThreeVector* theCenter = new G4ThreeVector(0,0,0);

  // The Pit.
  G4Material* theAir=new G4Material("N2-air", 14.0, 28.0 * g / mole,
				     0.0001 *g/cm3);
  G4Box* aBox = new G4Box("aBox",10*m, 10*m, 10*m);
  G4LogicalVolume*  theLogicalPit =
    new G4LogicalVolume(aBox, theAir,
			"thePit", NULL, NULL, NULL);
  theLogicalPit->SetVisAttributes (G4VisAttributes::Invisible);
  G4PVPlacement* thePit =
    new G4PVPlacement(    theNull,
			  *theCenter,
			  "thePit",
			  theLogicalPit,    
			  0, NULL, 0);

  // The Mother of Parametrised.
  G4VSolid * anotherBox = new G4Box("aBox",5*m, 5*m, 5*m);
  G4LogicalVolume*  theLogicalMother =
    new G4LogicalVolume(anotherBox, theAir,
			"aTest", NULL, NULL,
			NULL);
  G4PVPlacement* theMother =
    new G4PVPlacement(    0, G4ThreeVector (-4 * m, 0 , 0),
			  "theMother",
			  theLogicalMother,    
			  thePit, NULL, 0);

  // Parametrized volume.
  G4Material * theCopper =
    new G4Material("Quaterhard-Cu",25.0, 63.0 * g / mole,
		   6.8*g/cm3);
  G4Box * theBox = new G4Box("aBox",1,1,1);
  G4LogicalVolume*  theLogicalBox =
    new G4LogicalVolume(theBox, theCopper,
			"aTest", NULL, NULL,
			NULL);
  G4VisAttributes * logicalPAttributes;
  logicalPAttributes = new G4VisAttributes(G4Colour(1.,0.,0.));
  theLogicalBox->SetVisAttributes(logicalPAttributes);
  ParametrisedBox* thePar = new ParametrisedBox;
  G4PVParameterised * theParametrized=
    new G4PVParameterised("theParametrizedtest",
                theLogicalBox,
                theMother,
                kYAxis,
                2,
                thePar);

  // The Mother of Replica.
  G4VSolid * box3 = new G4Box("aBox", 1 * m, 2 * m, 3 * m);
  G4LogicalVolume*  theLogicalRMother =
    new G4LogicalVolume(box3, theAir,
			"aTest", NULL, NULL,
			NULL);
  G4PVPlacement* theRMother =
    new G4PVPlacement(    0, G4ThreeVector (4 * m, 0, 0),
			  "theRMother",
			  theLogicalRMother,    
			  thePit, NULL, 0);

  // Replicated volume.
  G4VSolid * RBox = new G4Box("aBox", 0.2 * m, 2 * m, 3 * m);
  G4LogicalVolume*  theLogicalRBox =
    new G4LogicalVolume(RBox, theCopper,
			"aTest", NULL, NULL,
			NULL);
  G4VisAttributes * logicalRAttributes;
  logicalRAttributes = new G4VisAttributes(G4Colour(0.,1.,0.));
  theLogicalRBox->SetVisAttributes(logicalRAttributes);
  G4VPhysicalVolume * theReplica=
    new G4PVReplica("theReplicatest",
                theLogicalRBox,
                theRMother,
                kXAxis,
                5, 0.4 * m, 0);

  // The Mother of Tubs Replica.
  G4VSolid * tubs = new G4Tubs ("aTubsBox", 0.5 * m, 1 * m, 2 * m,
				90. * deg, 180. * deg);
  G4LogicalVolume*  theLogicalTMother =
    new G4LogicalVolume(tubs, theAir,
			"aTest", NULL, NULL,
			NULL);
  G4PVPlacement* theTMother =
    new G4PVPlacement(    0, G4ThreeVector (4 * m, 4 * m, 0),
			  "theTMother",
			  theLogicalTMother,    
			  thePit, NULL, 0);

  // Replicated volume.
  G4VSolid * tubsPart = new G4Tubs ("aTubsPart", 0.5 * m, 1 * m, 2 * m,
				    0., 36. * deg);
  G4LogicalVolume*  theLogicalTubsPart =
    new G4LogicalVolume(tubsPart, theCopper,
			"aTest", NULL, NULL,
			NULL);
  G4VisAttributes * logicalTAttributes;
  logicalTAttributes = new G4VisAttributes(G4Colour(0.,0.,1.));
  theLogicalTubsPart->SetVisAttributes(logicalTAttributes);
  G4VPhysicalVolume * theTReplica=
    new G4PVReplica("theTReplicatest",
                theLogicalTubsPart,
                theTMother,
                kPhi,
                5, 36. * deg, 90. * deg);

  /*************** Bloggs's box to workaround bug in G4SmartVoxelHeader.

  G4Box * thebloggsBox = new G4Box("aBox",5*m, 4*m, 2*m);
  G4LogicalVolume*  theLogicalbloggs =
    new G4LogicalVolume(thebloggsBox, theCopper,
			"aTest", NULL, NULL,
			NULL);
  G4VisAttributes * bloggsAttributes = 
    new G4VisAttributes(G4Colour(0.,1.,0.));
  theLogicalbloggs->SetVisAttributes(bloggsAttributes);
  G4PVPlacement* bloggs =
    new G4PVPlacement(    theNull,
			  G4ThreeVector (10*m, -30*m, -40*m),
			  "bloggs",
			  theLogicalbloggs,
			  thePitPosition,    
			  NULL, 0);
			  ***************/

  return thePit;
}
