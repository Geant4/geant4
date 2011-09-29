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
// $Id: Tst202DetectorConstruction.cc,v 1.6 2010-11-30 17:20:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#include "Tst202DetectorConstruction.hh"

#include "Tst202DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4Orb.hh"
#include "G4Trap.hh"
#include "G4EllipticalTube.hh"
#include "G4EllipticalCone.hh"
#include "G4Polyhedra.hh"
#include "G4Polycone.hh"
#include "G4Tet.hh"
#include "G4TwistedTubs.hh"
#include "G4TwistedBox.hh"
#include "G4TwistedTrd.hh"
#include "G4GenericTrap.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4QuadrangularFacet.hh"
#include "G4TwistedTrap.hh"
#include "G4Paraboloid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4DisplacedSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4TwoVector.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4VisExtent.hh"
#include "G4RunManager.hh"
#include "G4Polyhedron.hh"

#include <vector>

Tst202DetectorConstruction::Tst202DetectorConstruction():
  verbosity(0)
{
  new Tst202DetectorMessenger(this);

  volumeSelection["calorimeter_boxes"]               = true;
  volumeSelection["box_in_box"]                      = true;
  volumeSelection["tracker_tube"]                    = true;
  volumeSelection["displaced_box"]                   = true;
  volumeSelection["Boolean_solids"]                  = true;
  volumeSelection["Tubes_etc"]                       = true;
  volumeSelection["Extra_placements"]                = true;
  volumeSelection["Sphere"]                          = true;
  volumeSelection["Polyhedra"]                       = true;
  volumeSelection["Polycone"]                        = true;
  volumeSelection["Orb"]                             = true;
  volumeSelection["Trapezoid"]                       = true;
  volumeSelection["Elliptical_Tube"]                 = true;
  volumeSelection["Elliptical_Cone"]                 = true;
  volumeSelection["G4Cons"]                          = true;
  volumeSelection["Simple_shared"]                   = true;
  volumeSelection["Boolean_for_logical_volume_test"] = true;
  volumeSelection["Radially_replicated_tube_sector"] = true;
  volumeSelection["tetrahedron"]                     = true;
  volumeSelection["Twisted_Tubs"]                    = true;
  volumeSelection["Twisted_Box"]                     = true;
  volumeSelection["Twisted_Trd"]                     = true;
  volumeSelection["Another_Twisted_Trd"]             = true;

  G4cout <<
    "Selectable volume items (with default value):"
    "\n  Change with /test202/select <name> true|false"
	 << G4endl;
  std::map<G4String,G4bool>::iterator i;
  for (i = volumeSelection.begin(); i != volumeSelection.end(); ++i) {
    G4cout << i->first
	   << " (" << std::boolalpha << i->second << std::noboolalpha << ")"
	   << G4endl;
  }

  expHall_x = 20. * m;
  expHall_y = 20. * m;
  expHall_z = 20. * m;

  calBox_x = 100.*cm;
  calBox_y = 50.*cm;
  calBox_z = 50.*cm;
  rotAngle = 30.*deg;
  calPos = 200.*cm;
  calMaterialName = "Pb";

  trackerRadius = 50.*cm;
  trackerHight = 100.*cm;
  trackerPos = -200.*cm;
}

Tst202DetectorConstruction::~Tst202DetectorConstruction()
{
  for (size_t i = 0; i < materialPointerStore.size(); i++) {
    delete materialPointerStore[i];
  }
}

void Tst202DetectorConstruction::SetVolume(G4String name, G4bool select) {
  if (volumeSelection.find(name) != volumeSelection.end()) {
    volumeSelection[name] = select;
    G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
  } else {
    if (name != "none") {
      G4cout << name << " not found.  Possibilities are:" << G4endl;
    }
  }
  std::map<G4String,G4bool>::iterator i;
  for (i = volumeSelection.begin(); i != volumeSelection.end(); ++i) {
    G4cout << i->first
	   << " (" << std::boolalpha << i->second << std::noboolalpha << ")"
	   << G4endl;
  }
}

G4VPhysicalVolume* Tst202DetectorConstruction::Construct()
{

  //------------------------------------------------------ materials

  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;

  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxygen", symbol="O", iz=8., a);

  density = 1.29e-03*g/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2);
  materialPointerStore.push_back(Air);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);

  a = 207.19*g/mole;
  density = 11.35*g/cm3;
  G4Material* Pb = new G4Material(name="Lead", z=82., a, density);
  materialPointerStore.push_back(Pb);

  a = 39.95*g/mole;
  density = 1.782e-03*g/cm3;
  G4Material* Ar = new G4Material(name="ArgonGas", z=18., a, density);
  materialPointerStore.push_back(Ar);

  a = 26.98*g/mole;
  density = 2.7*g/cm3;
  G4Material* Al = new G4Material(name="Aluminum", z=13., a, density);
  materialPointerStore.push_back(Al);

  a = 55.85*g/mole;
  density = 7.87*g/cm3;
  G4Material* Fe = new G4Material(name="Iron", z=26., a, density);
  materialPointerStore.push_back(Fe);

  //------------------------------------------------------ volumes

  //------------------------------ experimental hall
  G4Box * experimentalHall_box
    = new G4Box("expHall_b",expHall_x,expHall_y,expHall_z);
  G4LogicalVolume * experimentalHall_log
    = new G4LogicalVolume(experimentalHall_box,Air,"expHall_L",0,0,0);
  //  G4VisAttributes * experimentalHallVisAtt
  //      = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  //  experimentalHallVisAtt->SetForceWireframe(true);
  //  experimentalHall_log->SetVisAttributes(experimentalHallVisAtt);
  experimentalHall_log -> SetVisAttributes (&G4VisAttributes::GetInvisible());
  G4VPhysicalVolume * experimentalHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),"expHall_P",
                        experimentalHall_log,0,false,0);

  //------------------------------ calorimeter boxes
  G4VPhysicalVolume* calo_phys = 0;
  G4LogicalVolume * calorimeter_log = 0;
  if (volumeSelection["calorimeter_boxes"]) {
    G4Box * calorimeter_box
      = new G4Box("calorimeter_b",calBox_x,calBox_y,calBox_z);
    G4Material* calMat;
    if(calMaterialName=="Pb")
      { calMat = Pb; }
    else if(calMaterialName=="Al")
      { calMat = Al; }
    else if(calMaterialName=="Fe")
      { calMat = Fe; }
    else
      { calMat = Air; }
    calorimeter_log
      = new G4LogicalVolume(calorimeter_box,calMat,"calo_L",0,0,0);
    G4VisAttributes * calorimeterVisAtt
      = new G4VisAttributes(G4Colour(0.,0.,1.,0.5));
    //  calorimeterVisAtt->SetForceWireframe(true);
    calorimeter_log->SetVisAttributes(calorimeterVisAtt);
    for(G4int i=0;i<3;i++)
      {
	G4RotationMatrix rm;
	rm.rotateZ(i*rotAngle);
	if (verbosity) rm.print(G4cout);
	calo_phys =
	  new G4PVPlacement
	  (G4Transform3D(rm,G4ThreeVector(0.*cm,i*calPos,0.*cm)),
	   "calo_phys",calorimeter_log,experimentalHall_phys,
	   false,i);
	if (verbosity)
	  calo_phys->GetObjectRotationValue().print(G4cout);
      }
  }

  //------------------------------ box in box
  if (volumeSelection["box_in_box"]) {
    // This is a test of the BooleanProcessor for the case of a
    // "hollow" subtraction, i.e., making a hollow box using
    // G4SubtractionSolid.
    G4String name = "BoxInBox";
    G4double halfOuter = 50.*cm;
    G4double thickness = 10.*cm;
    G4double halfInner = halfOuter - thickness;
    G4VSolid* outer = new G4Box(name+"Outer", halfOuter, halfOuter, halfOuter);
    G4VSolid* inner = new G4Box(name+"Inner", halfInner, halfInner, halfInner);
    G4VSolid* hollow = new G4SubtractionSolid(name, outer, inner, G4Transform3D());
    G4LogicalVolume* log = new G4LogicalVolume(hollow, Al, name);
    new G4PVPlacement(G4Translate3D(0., 600.*cm, 0.),
		      log, name, experimentalHall_log, false, 0);
  }

  //------------------------------ tracker tube
  if (volumeSelection["tracker_tube"]) {
    G4Tubs * tracker_tube
      = new G4Tubs("tracker_tube",0.*cm,trackerRadius,trackerHight,
		   0.*deg,360.*deg);
    /*
      G4cout << "Tracker tube polyhedron:\n"
      << (HepPolyhedron)(*(tracker_tube->GetPolyhedron())) << G4endl;
    */
    G4LogicalVolume * tracker_log
      = new G4LogicalVolume(tracker_tube,Ar,"tracker_L",0,0,0);
    G4VisAttributes * trackerVisAtt
      = new G4VisAttributes(G4Colour(0.,0.,1.));
    //  trackerVisAtt->SetForceWireframe(true);
    tracker_log->SetVisAttributes(trackerVisAtt);
    //////////////............
    G4RotationMatrix* tracker_rm = new G4RotationMatrix;
    tracker_rm->rotateY(-30.*deg);
    if (verbosity)
      tracker_rm->print(G4cout);
    G4VPhysicalVolume* tracker_phys =
      //new G4PVPlacement(tracker_rm,G4ThreeVector(0.*cm,trackerPos,200.*cm),
      new G4PVPlacement
      (G4Transform3D(*tracker_rm,G4ThreeVector(0.*cm,trackerPos,200.*cm)),
       //////////////............
       //new G4PVPlacement(0,G4ThreeVector(0.*cm,trackerPos,0.*cm),
       "tracker_phys",tracker_log,experimentalHall_phys,
       false,0);
    if (verbosity)
      tracker_phys->GetObjectRotationValue().print(G4cout);
  }

  //------------------------------ displaced box
  if (volumeSelection["displaced_box"]) {
    G4Box * undisplaced_box =
      new G4Box("undisplaced_box",30.*cm,50.*cm,70.*cm);
    G4RotationMatrix rm_db;
    rm_db.rotateZ(20.*deg);
    G4DisplacedSolid* displaced_box = new G4DisplacedSolid
      ("displaced_box",undisplaced_box,
       G4Transform3D(rm_db,G4ThreeVector(200.*cm,0.,0.)));
    if (verbosity)
      G4cout << "Displaced box extent:\n" << displaced_box->GetExtent()
	     << G4endl;
    G4LogicalVolume * displaced_box_log
      = new G4LogicalVolume(displaced_box,Ar,"displaced_box_L",0,0,0);
    new G4PVPlacement(0,G4ThreeVector(0.*cm,-200.*cm,0.*cm),
		      "displaced_box_phys",displaced_box_log,
		      experimentalHall_phys,
		      false,0);
  }

  //-------------------------------------------- Boolean solids

  if (volumeSelection["Boolean_solids"]) {
    G4Tubs* cylinder1 = new G4Tubs("Cylinder #1",20*cm,50*cm,30*cm,0,2*pi);
    G4Box* box1 = new G4Box("Box #1",20*cm,30*cm,40*cm);
    G4Box* box2 = new G4Box("Box #2",10*cm,20*cm,35*cm);
    G4RotationMatrix* rm1 = new G4RotationMatrix;
    rm1->rotateZ(20*deg);
    G4RotationMatrix* rm2 = new G4RotationMatrix;
    rm2->rotateZ(60*deg);

    G4IntersectionSolid* cyl1Ibox1 =
      new G4IntersectionSolid("cylinder1-intersection-box1", cylinder1, box1,
			      rm1,G4ThreeVector(30.01*cm,30.01*cm,0.01*cm));
    if (verbosity)
      G4cout << "cylinder1-intersection-box1 extent:\n"
	     << cyl1Ibox1->GetExtent() << G4endl;
    G4IntersectionSolid* cyl1Ibox1Ibox2 =
      new G4IntersectionSolid
      ("cylinder1-intersection-box1-intersection-box2", cyl1Ibox1, box2,
       rm2,G4ThreeVector(0.01,40.01*cm,0.01*cm));
    if (verbosity)
      G4cout << "cylinder1-intersection-box1-intersection-box2 extent:\n"
	     << cyl1Ibox1Ibox2->GetExtent() << G4endl;
    G4LogicalVolume * intersection_log
      = new G4LogicalVolume(cyl1Ibox1Ibox2,Ar,"intersection_L",0,0,0);
    const G4VisAttributes* bool_red =
      new G4VisAttributes(G4Colour(1.,0.,0.));
    intersection_log->SetVisAttributes(bool_red);
    new G4PVPlacement
      (0,G4ThreeVector(100.*cm,-50*cm,0.*cm),
       "intersection_phys",intersection_log,experimentalHall_phys,
       true,0);

    G4SubtractionSolid* cyl1Sbox1 =
      new G4SubtractionSolid("cylinder1-subtraction-box1", cylinder1, box1,
			     rm1,G4ThreeVector(30.*cm,30.*cm,1.*cm));
    G4SubtractionSolid* cyl1Sbox1Sbox2 =
      new G4SubtractionSolid
      ("cylinder1-subtraction-box1-subtraction-box2", cyl1Sbox1, box2,
       rm2,G4ThreeVector(0.,40*cm,2.*cm));
    G4LogicalVolume * subtraction_log
      = new G4LogicalVolume(cyl1Sbox1Sbox2,Ar,"subtraction_L",0,0,0);
    const G4VisAttributes* bool_green =
      new G4VisAttributes(G4Colour(0.,1.,0.));
    subtraction_log->SetVisAttributes(bool_green);
    new G4PVPlacement
      (0,G4ThreeVector(200.*cm,-50*cm,0.*cm),
       "subtraction_phys",subtraction_log,experimentalHall_phys,
       true,0);

    G4UnionSolid* cyl1Ubox1 =
      new G4UnionSolid("cylinder1-union-box1", cylinder1, box1,
		       rm1,G4ThreeVector(30.*cm,30.*cm,1.*cm));
    G4UnionSolid* cyl1Ubox1Ubox2 =
      new G4UnionSolid("cylinder1-union-box1-union-box2", cyl1Ubox1, box2,
		       rm2,G4ThreeVector(0.,40*cm,2.*cm));
    G4LogicalVolume * union_log
      = new G4LogicalVolume(cyl1Ubox1Ubox2,Ar,"union_L",0,0,0);
    const G4VisAttributes* bool_blue =
      new G4VisAttributes(G4Colour(0.,0.,1.));
    union_log->SetVisAttributes(bool_blue);
    new G4PVPlacement(0,G4ThreeVector(350.*cm,-50*cm,0.*cm),
		      "union_phys",union_log,experimentalHall_phys,
		      true,0);
  }

  //----------- Tubes, replicas(!?) and daughter boxes

  G4LogicalVolume * divided_tube_inset_log = 0;
  G4LogicalVolume * grand_daughter_box2_log = 0;
  if (volumeSelection["Tubes_etc"]) {
    const G4double eps = 10 * mm;
    const G4double alp = 10 * mrad;

    G4VisAttributes * grey
      = new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.1));

    G4double tube_dPhi = twopi;
    G4Tubs* tube = new G4Tubs("tube",20*cm,50*cm,30*cm,0.,tube_dPhi);
    G4LogicalVolume * tube_log
      = new G4LogicalVolume(tube,Ar,"tube_L",0,0,0);
    //G4VisAttributes * tube_VisAtt
    //  = new G4VisAttributes(G4Colour(0.,1.,0.,0.1));
    //  tube_log->SetVisAttributes(tube_VisAtt);
    tube_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    new G4PVPlacement(0,G4ThreeVector(-200.*cm,0.,0.*cm),
		      "tube_phys",tube_log,experimentalHall_phys,
		      false,0);

    G4double divided_tube_dPhi = tube_dPhi / 6.;
    G4Tubs* divided_tube = new G4Tubs
      ("divided_tube",
       20*cm,50*cm,30*cm,-divided_tube_dPhi/2.,divided_tube_dPhi);
    G4LogicalVolume * divided_tube_log
      = new G4LogicalVolume(divided_tube,Ar,"divided_tube_L",0,0,0);
    divided_tube_log->SetVisAttributes(grey);
    new G4PVReplica("divided_tube_phys",divided_tube_log,tube_log,
		    kPhi,6,divided_tube_dPhi);
    /************ 
  G4int iCopy;
  for (iCopy = 0; iCopy < 6; iCopy++) {
    new G4PVPlacement
      (G4Transform3D
       (G4RotationMatrix().rotateZ(divided_tube_dPhi/2.+iCopy*pi/3.),
	G4ThreeVector()),
       divided_tube_log,"divided_tube_phys",tube_log,
       false,iCopy);
  }
    *********/

    G4double divided_tube_inset_dPhi = divided_tube_dPhi - 2. * alp;
    G4Tubs* divided_tube_inset = new G4Tubs
      ("divided_tube_inset",
       20*cm+eps,50*cm-eps,30*cm-eps,
       -divided_tube_inset_dPhi/2.,
       divided_tube_inset_dPhi);
    divided_tube_inset_log = new G4LogicalVolume
      (divided_tube_inset,Ar,"divided_tube_inset_L",0,0,0);
    //G4VisAttributes * divided_tube_inset_VisAtt
    //  = new G4VisAttributes(G4Colour(1.,0.,0.,0.2));
    //  divided_tube_inset_log->SetVisAttributes(divided_tube_inset_VisAtt);
    divided_tube_inset_log->SetVisAttributes(G4VisAttributes::Invisible);
    new G4PVPlacement(0,G4ThreeVector(),
		      divided_tube_inset_log,"divided_tube_inset_phys",
		      divided_tube_log,false,0);

    G4double sub_divided_tube_dPhi = divided_tube_inset_dPhi / 4.;
    G4Tubs* sub_divided_tube = new G4Tubs
      ("sub_divided_tube",
       20*cm+eps,50*cm-eps,30*cm-eps,
       -sub_divided_tube_dPhi/2.,sub_divided_tube_dPhi);
    G4LogicalVolume * sub_divided_tube_log
      = new G4LogicalVolume(sub_divided_tube,Ar,"sub_divided_tube_L",0,0,0);
    sub_divided_tube_log->SetVisAttributes(grey);
    new G4PVReplica("sub_divided_tube_phys",
		    sub_divided_tube_log,divided_tube_inset_log,
		    kPhi,4,sub_divided_tube_dPhi,-divided_tube_inset_dPhi/2.);
    /************ 
  for (iCopy = 0; iCopy < 4; iCopy++) {
    new G4PVPlacement
      (G4Transform3D
       (G4RotationMatrix().rotateZ
	(-divided_tube_inset_dPhi/2.
	 +(iCopy+0.5)*sub_divided_tube_dPhi),
        G4ThreeVector()),
       sub_divided_tube_log,"sub_divided_tube_phys",divided_tube_inset_log,
       false,iCopy);
  }
    ************/
  
    G4Box* daughter_box = new G4Box("daughter_box",4.*cm,3.*cm,25.*cm);
    G4LogicalVolume * daughter_box_log
      = new G4LogicalVolume(daughter_box,Ar,"daughter_box_L",0,0,0);
    G4VisAttributes * daughter_box_VisAtt
      = new G4VisAttributes(G4Colour(0.,0.,1.,0.3));
    daughter_box_log->SetVisAttributes(daughter_box_VisAtt);
    G4Box* grand_daughter_box = new G4Box("grand_daughter_box",1*cm,2*cm,5*cm);
    G4LogicalVolume * grand_daughter_box_log
      = new G4LogicalVolume
      (grand_daughter_box,Ar,"grand_daughter_box_L",0,0,0);
    G4VisAttributes * grand_daughter_box_VisAtt
      = new G4VisAttributes(G4Colour(1.,1.,0.));
    grand_daughter_box_log->SetVisAttributes(grand_daughter_box_VisAtt);
    G4Box* grand_daughter_box2 =
      new G4Box("grand_daughter_box2",1*cm,2*cm,5*cm);
    grand_daughter_box2_log = new G4LogicalVolume
      (grand_daughter_box2,Ar,"grand_daughter_box2_L",0,0,0);
    G4VisAttributes * grand_daughter_box2_VisAtt
      = new G4VisAttributes(G4Colour(1.,0.,1.));
    grand_daughter_box2_log->SetVisAttributes(grand_daughter_box2_VisAtt);
    new G4PVPlacement(0,G4ThreeVector(-2*cm,0.,0.),
		      grand_daughter_box_log,"grand_daughter_box_phys",
		      daughter_box_log,false,0);
    new G4PVPlacement(0,G4ThreeVector(2*cm,0.,0.),
		      grand_daughter_box2_log,"grand_daughter_box2_phys",
		      daughter_box_log,false,0);
    new G4PVPlacement(0,G4ThreeVector(40*cm,0.,0.),
		      daughter_box_log,"daughter_box_phys",
		      sub_divided_tube_log,false,0);
  }

  //------------------------------------------------ Extra placements
  // For good measure, as a test of drawn volume path, place in one of
  // the earlier volumes...
  if (volumeSelection["Extra_placements"]) {
    new G4PVPlacement(0,G4ThreeVector(),
		      "divided_tube_inset_phys",divided_tube_inset_log,
		      calo_phys,false,0);  // Place in PV.
    new G4PVPlacement(0,G4ThreeVector(),
		      grand_daughter_box2_log,"grand_daughter_box2_phys",
		      calorimeter_log,false,0);  // Place in LV.  Same effect.
  }

  //-------------------------------------------- Sphere

  if (volumeSelection["Sphere"]) {
    G4Sphere* PD_vol_crystal
      = new G4Sphere("Test_Sphere",
                     100.*cm,           // inner radius
                     200.*cm,           // outer radius
                     0.,                // start phi
                     90.*degree,         // delta phi
                     0.,                // start theta
                     90.*degree          // delta theta
                     );

    G4LogicalVolume* PD_log_crystal
      = new G4LogicalVolume(PD_vol_crystal,Ar,"Test_Sphere");

    G4VisAttributes * PD_att_crystal
      = new G4VisAttributes(G4Colour(1.,0.,1.));
    PD_att_crystal->SetForceAuxEdgeVisible(true);
    PD_att_crystal->SetForceLineSegmentsPerCircle(100);
    PD_log_crystal->SetVisAttributes(PD_att_crystal);

    G4RotationMatrix rm;

    new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(200.*cm,300.*cm,0)),
		      "PD_physical", PD_log_crystal,
		      experimentalHall_phys,false,0);
  }

//-------------------------------------------- Polyhedra and Polycone

  if (volumeSelection["Polyhedra"]) {
    const G4int numRZ = 10;
    G4double polyhedra_r[] = {0,5,3,4,9,9,3,3,2,0};
    G4double polyhedra_z[] = {0,1,2,3,0,5,4,3,2,1};
    for (int i = 0; i < numRZ; ++i) {
      polyhedra_r[i] *= 10*cm;
      polyhedra_z[i] *= 10.*cm;
    }

    G4Polyhedra* polyhedra_solid
      = new G4Polyhedra("Polyhedra_Test",
			0.,270.*deg,6,numRZ,polyhedra_r,polyhedra_z);
    //0.,twopi,6,numRZ,polyhedra_r,polyhedra_z);
    if (verbosity)
      G4cout << polyhedra_solid->StreamInfo(G4cout) << G4endl;

    G4LogicalVolume* polyhedra_log
      = new G4LogicalVolume(polyhedra_solid,Ar,"Polyhedra_Test");

    G4VisAttributes * polyhedra_atts
      = new G4VisAttributes(G4Colour(0.,1.,1.));
    polyhedra_log->SetVisAttributes(polyhedra_atts);

    G4RotationMatrix rm;

    new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(200.*cm,100.*cm,0)),
		      "Polyhedra_Test", polyhedra_log,
		      experimentalHall_phys,false,0);
  }

  if (volumeSelection["Polycone"]) {
    const G4int numRZ1 = 10;
    G4double polycone_r[] = {1,5,3,4,9,9,3,3,2,1};
    G4double polycone_z[] = {0,1,2,3,0,5,4,3,2,1};
    for (int i = 0; i < numRZ1; ++i) {
      polycone_r[i] *= 10*cm;
      polycone_z[i] *= 10.*cm;
    }

    G4Polycone* polycone_solid
      = new G4Polycone("Polycone_Test",
		       0.*deg,270.*deg,numRZ1,polycone_r,polycone_z);
    //0.,twopi,numRZ,polycone_r,polycone_z);
    if (verbosity)
      G4cout << polycone_solid->StreamInfo(G4cout) << G4endl;

    G4LogicalVolume* polycone_log
      = new G4LogicalVolume(polycone_solid,Ar,"Polycone_Test");

    G4VisAttributes * polycone_atts
      = new G4VisAttributes(G4Colour(1.,0.5,0.5));
    polycone_log->SetVisAttributes(polycone_atts);

    G4RotationMatrix rm;

    new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(400.*cm,120.*cm,0)),
		      "Polycone_Test", polycone_log,
		      experimentalHall_phys,false,0);
  }

  //-------------------------------------------- Orb

  if (volumeSelection["Orb"]) {
    G4Orb* orb 
      = new G4Orb("Test orb", 100.*cm);

    G4LogicalVolume* orb_log
      = new G4LogicalVolume(orb, Ar,"Test_orb");

    G4VisAttributes * orb_att
      = new G4VisAttributes(G4Colour(1.,0.,1.));
    orb_att->SetForceAuxEdgeVisible(true);
    orb_log->SetVisAttributes(orb_att);

    new G4PVPlacement(G4Translate3D(-300.*cm,200.*cm,0),
		      "Test orb", orb_log, 
		      experimentalHall_phys,false,0);
  }

  //-------------------------------------------- Trapezoid

  if (volumeSelection["Trapezoid"]) {
    G4Trap* trap1_solid = new G4Trap
      ("trap1_solid",
       50.*cm,  //      pDz     Half-length along the z-axis
       0.,      //      pTheta  Polar angle of the line joining the centres of the faces
       //              at -/+pDz
       0.,      //      pPhi    Azimuthal angle of the line joing the centre of the face at
       //              -pDz to the centre of the face at +pDz
       40.*cm,  //      pDy1     Half-length along y of the face at -pDz
       30.*cm,  //      pDx1    Half-length along x of the side at y=-pDy1 of the face at -pDz
       30.*cm,  //      pDx2    Half-length along x of the side at y=+pDy1 of the face at -pDz
       20.*deg, //      pAlp1   Angle with respect to the y axis from the centre of the side
       //              at y=-pDy1 to the centre at y=+pDy1 of the face at -pDz
       40.*cm,  //      pDy2     Half-length along y of the face at +pDz
       30.*cm,  //      pDx3    Half-length along x of the side at y=-pDy2 of the face at +pDz
       30.*cm,  //      pDx4    Half-length along x of the side at y=+pDy2 of the face at +pDz
       20.*deg  //      pAlp2   Angle with respect to the y axis from the centre of the side
       //              at y=-pDy2 to the centre at y=+pDy2 of the face at +pDz
       );
    G4LogicalVolume* trap1_log = new G4LogicalVolume (trap1_solid,Ar,"trap1_log");

    G4RotationMatrix rm;

    new G4PVPlacement
      (G4Transform3D(rm,G4ThreeVector(-200.*cm,200.*cm,-200.*cm)),
       "trap1_phys", trap1_log,
       experimentalHall_phys,false,0);

    G4Trap* trap2_solid = new G4Trap
      ("trap2_solid",
       50.*cm,  //      pDz     Half-length along the z-axis
       20.*deg, //      pTheta  Polar angle of the line joining the centres of the faces
       //              at -/+pDz
       90.*deg, //      pPhi    Azimuthal angle of the line joing the centre of the face at
       //              -pDz to the centre of the face at +pDz
       40.*cm,  //      pDy1     Half-length along y of the face at -pDz
       30.*cm,  //      pDx1    Half-length along x of the side at y=-pDy1 of the face at -pDz
       30.*cm,  //      pDx2    Half-length along x of the side at y=+pDy1 of the face at -pDz
       0.,      //      pAlp1   Angle with respect to the y axis from the centre of the side
       //              at y=-pDy1 to the centre at y=+pDy1 of the face at -pDz
       40.*cm,  //      pDy2     Half-length along y of the face at +pDz
       30.*cm,  //      pDx3    Half-length along x of the side at y=-pDy2 of the face at +pDz
       30.*cm,  //      pDx4    Half-length along x of the side at y=+pDy2 of the face at +pDz
       0.       //      pAlp2   Angle with respect to the y axis from the centre of the side
       //              at y=-pDy2 to the centre at y=+pDy2 of the face at +pDz
       );
    G4LogicalVolume* trap2_log = new G4LogicalVolume (trap2_solid,Ar,"trap2_log");
    new G4PVPlacement
      (G4Transform3D(rm,G4ThreeVector(-200.*cm,400.*cm,-200.*cm)),
       "trap2_phys", trap2_log,
       experimentalHall_phys,false,0);
  }

  //-------------------------------------------- Elliptical Tube
  if (volumeSelection["Elliptical_Tube"]) {
    G4VSolid* eTube = new G4EllipticalTube("e-tube",100.*cm,50.*cm,100.*cm);
    G4LogicalVolume* eTubeLog = new G4LogicalVolume(eTube,Ar,"e-tube-log");
    new G4PVPlacement(G4Translate3D(G4ThreeVector(-400.*cm,0.,0)),
		      "e-tube-phys", eTubeLog,
		      experimentalHall_phys,false,0);
  }

  //-------------------------------------------- Elliptical Cone
  if (volumeSelection["Elliptical_Cone"]) {
    //Creating EllipticalCone with Dx=50*cm and Dy=100.*cm at z=-50.*cm
    
    G4VSolid* eCone = new G4EllipticalCone("e-cone",1./3.,2./3.,100.*cm,50.*cm);
    //G4VSolid* eCone = new G4EllipticalCone("e-cone",50.*cm,100.*cm,100.*cm,50.*cm);
    //G4VSolid* eCone = new G4EllipticalCone("e-cone",1.*mm,0.5*mm,40.*mm,20.*mm);
    G4LogicalVolume* eConeLog = new G4LogicalVolume(eCone,Ar,"e-cone-log");
    new G4PVPlacement(G4Translate3D(G4ThreeVector(-500.*cm,400.*cm,0)),
		      "e-cone-phys", eConeLog,
		      experimentalHall_phys,false,0);
  }

  //-------------------------------------------- G4Cons
  if (volumeSelection["G4Cons"]) {
    G4VSolid* eCons = new G4Cons
      ("e-cons",50.*cm,70.*cm,100.*cm,140.*cm,200.*cm,0,twopi);
    G4LogicalVolume* eConsLog = new G4LogicalVolume(eCons,Ar,"e-cons-log");
    new G4PVPlacement(G4Translate3D(G4ThreeVector(-300.*cm,-200.*cm,60*cm)),
		      "e-cons-phys", eConsLog,
		      experimentalHall_phys,false,0);
  }

  //--------------------------- Simple shared logical volume tree

  if (volumeSelection["Simple_shared"]) {
    G4VSolid* boxTwo = new G4Box("Box2", 100.*cm, 100.*cm, 100.*cm);
    G4LogicalVolume* log2 = new G4LogicalVolume(boxTwo,Ar,"Log2");
    new G4PVPlacement(G4Translate3D(G4ThreeVector(600.*cm,0.,0.)),
		      "B",log2,
		      experimentalHall_phys,false,0);
    new G4PVPlacement(G4Translate3D(G4ThreeVector(900.*cm,0.,0.)),
		      "C",log2,
		      experimentalHall_phys,false,0);
    G4VSolid* boxThree = new G4Box("Box3", 50.*cm, 50.*cm, 50.*cm);
    G4LogicalVolume* log3 = new G4LogicalVolume(boxThree,Ar,"Log3");
    new G4PVPlacement(G4Translate3D(),
		      log3,"D",
		      log2,false,0);
    G4VSolid* boxFour = new G4Box("Box3", 30.*cm, 30.*cm, 30.*cm);
    G4LogicalVolume* log4 = new G4LogicalVolume(boxFour,Ar,"Log4");
    new G4PVPlacement(G4Translate3D(),
		      log4,"E",
		      log3,false,0);
  }

  //--------------------------- Boolean for logical volume test

  if (volumeSelection["Boolean_for_logical_volume_test"]) {
    G4double alSLayer_x = 30.*cm;
    G4double alSLayer_y = 50.*cm;
    G4double alSLayer_z = 2.5*cm;

    G4Box* box1B = new
      G4Box("alSLayer_box",alSLayer_x/2.,alSLayer_y/2.,alSLayer_z/2.);

    G4Tubs* Cylinder1 = new G4Tubs
      ("Cylinder#1",0.*mm,alSLayer_x/2.,alSLayer_x/2.,0.,2.*CLHEP::pi); 

    G4UnionSolid* b1UnionC1 = new G4UnionSolid("Box+Cylinder", box1B,
					       Cylinder1); 

    G4LogicalVolume* alSpaceCraft_log = new
      G4LogicalVolume(b1UnionC1,Al,"alLayer_log",0,0,0);  

    new G4PVPlacement(0,G4ThreeVector(900.*cm, 200.*cm, 0.),alSpaceCraft_log,"alLayer_phys",experimentalHall_log,false,0);
    // new G4PVPlacement(G4Translate3D(G4ThreeVector(900.*cm, 200.*cm, 0.)),"alLayer_phys",alSpaceCraft_log,experimentalHall_log,false,0);
  }

  //----------- Radially replicated tube sector

  if (volumeSelection["Radially_replicated_tube_sector"]) {
    G4double rMin = 50.*cm;
    G4double DeltaR = 50.*cm;
    G4VSolid* rrTubs = new G4Tubs
      ("rrTubs",rMin,rMin + DeltaR,200*cm,180*deg,90*deg);
    G4LogicalVolume* rrTubsLog = new G4LogicalVolume
      (rrTubs,Ar,"rrTubs-log");
    rrTubsLog->SetVisAttributes(G4VisAttributes::Invisible);
    new G4PVPlacement(G4Translate3D(G4ThreeVector(400.*cm,-200.*cm,0)),
		      "rrTubs-phys", rrTubsLog,
		      experimentalHall_phys,false,0);
    G4double deltaR = DeltaR / 6.;
    G4Tubs* drrTubs = new G4Tubs
      ("drrTubs",rMin,rMin + deltaR,200*cm,180*deg,90*deg);
    G4LogicalVolume * drrTubsLog = new G4LogicalVolume
      (drrTubs,Ar,"drrTubs-log");
    drrTubsLog->SetVisAttributes(G4Colour::Red());
    new G4PVReplica("drrTubs-phys",drrTubsLog,rrTubsLog,
		    kRho,6,deltaR,rMin);
  }


  //----------- tetrahedron

  if (volumeSelection["tetrahedron"]) {
    G4VSolid* tet = new G4Tet
      ("tet",
       G4ThreeVector(),
       G4ThreeVector(0.,100.*cm,0.),
       G4ThreeVector(0.,0.,100.*cm),
       G4ThreeVector(100*cm,0.,0.));
    G4LogicalVolume* tet_log = new G4LogicalVolume
      (tet,Ar,"tet-log");
    tet_log->SetVisAttributes(G4VisAttributes(G4Colour(0.,1.,0.)));
    new G4PVPlacement
      (G4Translate3D(G4ThreeVector(300.*cm,-400.*cm,0.)),
       "tet-phys", tet_log,
       experimentalHall_phys,false,0);
  }

  //----------- Twisted solids

  G4VSolid* aSolid;
  G4LogicalVolume* aLog;
  G4double fTrackerR1 ;   // r1
  G4double fTrackerR2 ;   // r2
  G4double fTrackerpDz  ;  // Full length of Tracker (pDz)
  G4double fTrackerpDx1 ;  // twisted Trapezoid
  G4double fTrackerpDx2 ;
  G4double fTrackerpDy1 ;
  G4double fTrackerpDy2 ;
  G4double fTwistAngle ;
  G4double fPhi ;

  G4double myScale = 5.;

  if (volumeSelection["Twisted_Tubs"]) {
    fTwistAngle = 90*deg ;
    fTrackerpDz = myScale*20*cm ;
    fTrackerR1  = myScale*7*cm ;
    fTrackerR2  = myScale*10*cm ;
    fPhi        = 180*deg ;

    aSolid = new G4TwistedTubs
      ("aTwistedTubs", fTwistAngle,
       fTrackerR1, fTrackerR2, fTrackerpDz, fPhi ) ;
    //G4cout << "aTwistedTubs: volume: "
    //       << aSolid->GetCubicVolume() << G4endl;
    aLog = new G4LogicalVolume
      (aSolid,Ar,"aTwistedTubs-log");
    aLog->SetVisAttributes(G4VisAttributes(G4Colour(1.,1.,0.)));
    new G4PVPlacement
      (G4Translate3D(G4ThreeVector(200.*cm,-400.*cm,0.)),
       "aTwistedTubs-phys", aLog, experimentalHall_phys,false,0);
  }

  if (volumeSelection["Twisted_Box"]) {
    fTwistAngle = 50*deg ;
    fTrackerpDx1 = myScale*4*cm ;
    fTrackerpDy1 = myScale*6*cm ;
    fTrackerpDz  = myScale*15*cm ;

    aSolid = new G4TwistedBox
      ("aTwistedBox",fTwistAngle,fTrackerpDx1,fTrackerpDy1,fTrackerpDz) ;
    aLog = new G4LogicalVolume
      (aSolid,Ar,"aTwistedBox-log");
    aLog->SetVisAttributes(G4VisAttributes(G4Colour(1.,1.,0.)));
    new G4PVPlacement
      (G4Translate3D(G4ThreeVector(100.*cm,-400.*cm,0.)),
       "aTwistedBox-phys", aLog, experimentalHall_phys,false,0);
  }

  if (volumeSelection["Twisted_Trd"]) {
    fTrackerpDx1 = myScale*4*cm ;
    fTrackerpDx2 = myScale*7*cm ;
    fTrackerpDy1 = myScale*2*cm ;
    fTrackerpDy2 = myScale*4*cm ;
    fTrackerpDz  = myScale*15*cm ;
    fTwistAngle = 50*deg ;

    aSolid = new G4TwistedTrd
      ("aTwistedTrd",fTrackerpDx1,fTrackerpDx2,fTrackerpDy1,fTrackerpDy2,
       fTrackerpDz,fTwistAngle);
    aLog = new G4LogicalVolume
      (aSolid,Ar,"aTwistedTrd-log");
    aLog->SetVisAttributes(G4VisAttributes(G4Colour(1.,1.,0.)));
    new G4PVPlacement
      (G4Translate3D(G4ThreeVector(000.*cm,-400.*cm,0.)),
       "aTwistedTrd-phys", aLog, experimentalHall_phys,false,0);
  }

  if (volumeSelection["Another_Twisted_Trd"]) {
    fTrackerpDx1 = myScale*4*cm ;
    fTrackerpDx2 = myScale*7*cm ;
    fTrackerpDy1 = myScale*2*cm ;
    fTrackerpDy2 = myScale*4*cm ;
    fTrackerpDz  = myScale*15*cm ;
    fTwistAngle = 50*deg ;

    aSolid = new G4TwistedTrd
      ("anotherTwistedTrd",fTrackerpDx1,fTrackerpDx2,fTrackerpDy1,fTrackerpDy2,
       fTrackerpDz,fTwistAngle);
    aLog = new G4LogicalVolume
      (aSolid,Ar,"anotherTwistedTrd-log");
    aLog->SetVisAttributes(G4VisAttributes(G4Colour(1.,1.,0.)));
    new G4PVPlacement
      (G4Translate3D(G4ThreeVector(-100.*cm,-400.*cm,0.)),
       "anotherTwistedTrd-phys", aLog, experimentalHall_phys,false,0);
  }

  aSolid = new G4Paraboloid("paraboloid",20.*cm,5.*cm,35.*cm);
  aLog = new G4LogicalVolume(aSolid,Ar,"paraboloid-log");
  new G4PVPlacement
    (G4Translate3D(G4ThreeVector(-200.*cm,-400.*cm,0.)),
       "paraboloid-phys", aLog, experimentalHall_phys,false,0);

  //--- Create a simple box using G4GenericTrap and make a Boolean operation
  std::vector<G4TwoVector> vertices2; 
  vertices2.push_back( G4TwoVector( 50.*cm, 50.*cm) );
  vertices2.push_back( G4TwoVector( 50.*cm,-50.*cm) );
  vertices2.push_back( G4TwoVector(-50.*cm,-50.*cm) );
  vertices2.push_back( G4TwoVector(-50.*cm, 50.*cm) );
  vertices2.push_back( G4TwoVector( 50.*cm, 50.*cm) );
  vertices2.push_back( G4TwoVector( 50.*cm,-50.*cm) );
  vertices2.push_back( G4TwoVector(-50.*cm,-50.*cm) );
  vertices2.push_back( G4TwoVector(-50.*cm, 50.*cm) );
  G4VSolid* testTrap = new G4GenericTrap( "TestTrap", 100.*cm, vertices2 );
  G4LogicalVolume* testTrap_log =
    new G4LogicalVolume( testTrap, Ar, "TestTrap-log" );
  new G4PVPlacement
    (G4Translate3D(G4ThreeVector(600.*cm,-300.*cm,0.)),
       "TestTrapCut-phys", testTrap_log, experimentalHall_phys, false, 0 );

  // Boolean subtraction
  G4VSolid* testCut = new G4Box( "TestCut", 10.*cm, 10.*cm, 10.*cm);
  G4VSolid* testTrapCut
    = new G4SubtractionSolid( "TestTrapCut", testTrap, testCut,
			      G4Translate3D(50.*cm, 50.*cm, 100.*cm) );
  G4LogicalVolume* testTrapCut_log =
    new G4LogicalVolume( testTrapCut, Ar, "TestTrapCut-log" );
  new G4PVPlacement
    (G4Translate3D(G4ThreeVector(800.*cm,-300.*cm,0.)),
       "TestTrapCut-phys", testTrapCut_log, experimentalHall_phys, false, 0 );

  //--- Create a simple tessellated solid
  G4double targetSize = 100*cm ;
  G4TriangularFacet *facet1 = new
  G4TriangularFacet (G4ThreeVector(-targetSize,-targetSize,        0.0),
                     G4ThreeVector(+targetSize,-targetSize,        0.0),
                     G4ThreeVector(        0.0,        0.0,+targetSize),
                     ABSOLUTE);
  G4TriangularFacet *facet2 = new
  G4TriangularFacet (G4ThreeVector(+targetSize,-targetSize,        0.0),
                     G4ThreeVector(+targetSize,+targetSize,        0.0),
                     G4ThreeVector(        0.0,        0.0,+targetSize),
                     ABSOLUTE);
  G4TriangularFacet *facet3 = new
  G4TriangularFacet (G4ThreeVector(+targetSize,+targetSize,        0.0),
                     G4ThreeVector(-targetSize,+targetSize,        0.0),
                     G4ThreeVector(        0.0,        0.0,+targetSize),
                     ABSOLUTE);
  G4TriangularFacet *facet4 = new
  G4TriangularFacet (G4ThreeVector(-targetSize,+targetSize,        0.0),
                     G4ThreeVector(-targetSize,-targetSize,        0.0),
                     G4ThreeVector(        0.0,        0.0,+targetSize),
                     ABSOLUTE);
  G4QuadrangularFacet *facet5 = new
  G4QuadrangularFacet (G4ThreeVector(-targetSize,-targetSize,      0.0),
                       G4ThreeVector(-targetSize,+targetSize,      0.0),
                       G4ThreeVector(+targetSize,+targetSize,      0.0),
                       G4ThreeVector(+targetSize,-targetSize,      0.0),
                       ABSOLUTE);
  G4TessellatedSolid* testTess = new G4TessellatedSolid("TestTess");
  testTess->AddFacet(facet1);
  testTess->AddFacet(facet2);
  testTess->AddFacet(facet3);
  testTess->AddFacet(facet4);
  testTess->AddFacet(facet5);
  testTess->SetSolidClosed(true);
  G4LogicalVolume* testTess_log =
    new G4LogicalVolume( testTess, Ar, "TestTess-log" );
  new G4PVPlacement
    (G4Translate3D(G4ThreeVector(1000.*cm,-300.*cm,0.)),
       "TestTessCut-phys", testTess_log, experimentalHall_phys, false, 0 );

  // Boolean subtraction
  G4VSolid* testCut1 = new G4Box( "TestCut", 100.*cm, 100.*cm, 100.*cm);
  G4VSolid* testTessCut
    = new G4SubtractionSolid( "TestTessCut", testTess, testCut1,
			      G4Translate3D(0., 0., 150.*cm) );
  G4LogicalVolume* testTessCut_log =
    new G4LogicalVolume( testTessCut, Ar, "TestTessCut-log" );
  new G4PVPlacement
    (G4Translate3D(G4ThreeVector(1300.*cm,-300.*cm,0.)),
       "TestTessCut-phys", testTessCut_log, experimentalHall_phys, false, 0 );

  //-------------------------------------------- return
  return experimentalHall_phys;
}
