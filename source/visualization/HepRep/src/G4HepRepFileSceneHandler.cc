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
// $Id: G4HepRepFileSceneHandler.cc 101714 2016-11-22 08:53:13Z gcosmo $
//
//
// Joseph Perl  27th January 2002
// A base class for a scene handler to export geometry and trajectories
// to the HepRep xml file format.

#include "G4HepRepFileSceneHandler.hh"
#include "G4HepRepFile.hh"
#include "G4HepRepMessenger.hh"
#include "G4UIcommand.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Version.hh"
#include "G4VSolid.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include "G4Polymarker.hh"
#include "G4Polyline.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Polyhedron.hh"
#include "G4VTrajectory.hh"
#include "G4VTrajectoryPoint.hh"
#include "G4TrajectoriesModel.hh"
#include "G4VHit.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4AttCheck.hh"
#include "G4VisManager.hh"
#include "G4VisTrajContext.hh"
#include "G4VTrajectoryModel.hh"

//HepRep
#include "G4HepRepFileXMLWriter.hh"

G4int G4HepRepFileSceneHandler::fSceneIdCount = 0;
// Counter for HepRep scene handlers.

G4HepRepFileSceneHandler::G4HepRepFileSceneHandler(G4VGraphicsSystem& system,
												   const G4String& name):
G4VSceneHandler(system, fSceneIdCount++, name)
{
	hepRepXMLWriter = ((G4HepRepFile*)(&system))->GetHepRepXMLWriter();
	fileCounter = 0;
	
	inPrimitives2D = false;
	warnedAbout3DText = false;
	warnedAbout2DMarkers = false;
	haveVisible = false;
	drawingTraj = false;
	doneInitTraj = false;
	drawingHit = false;
	doneInitHit = false;
        trajContext = 0;
        trajAttValues = 0;
        trajAttDefs = 0;
        hitAttValues = 0;
        hitAttDefs = 0;
}


G4HepRepFileSceneHandler::~G4HepRepFileSceneHandler() {}


void G4HepRepFileSceneHandler::BeginModeling() {
    G4VisManager* visManager = G4VisManager::GetInstance();
    const G4VTrajectoryModel* model = visManager->CurrentTrajDrawModel();
    trajContext = & model->GetContext();
	
	G4VSceneHandler::BeginModeling();  // Required: see G4VSceneHandler.hh.
}


void G4HepRepFileSceneHandler::EndModeling() {
	G4VSceneHandler::EndModeling();  // Required: see G4VSceneHandler.hh.
}

void G4HepRepFileSceneHandler::BeginPrimitives2D
(const G4Transform3D& objectTransformation) {
#ifdef G4HEPREPFILEDEBUG
	G4cout << "G4HepRepFileSceneHandler::BeginPrimitives2D() " << G4endl;
#endif
	inPrimitives2D = true;
	G4VSceneHandler::BeginPrimitives2D(objectTransformation);
}

void G4HepRepFileSceneHandler::EndPrimitives2D () {
#ifdef G4HEPREPFILEDEBUG
	G4cout << "G4HepRepFileSceneHandler::EndPrimitives2D() " << G4endl;
#endif
	G4VSceneHandler::EndPrimitives2D();
	inPrimitives2D = false;
}


#ifdef G4HEPREPFILEDEBUG
void G4HepRepFileSceneHandler::PrintThings() {
	G4cout <<
    "  with transformation "
	<< fObjectTransformation;
	G4PhysicalVolumeModel* pPVModel =
		dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
	if (pPVModel) {
		G4VPhysicalVolume* pCurrentPV = pPVModel->GetCurrentPV();
		G4LogicalVolume* pCurrentLV = pPVModel->GetCurrentLV();
		G4int currentDepth = pPVModel->GetCurrentDepth();
		G4cout <<
			"\n  current physical volume: "
			<< pCurrentPV->GetName() <<
			"\n  current logical volume: "
			<< pCurrentLV->GetName() <<
			"\n  current depth of geometry tree: "
			<< currentDepth;
	}
	G4cout << G4endl;
}
#endif


void G4HepRepFileSceneHandler::AddSolid(const G4Box& box) {
#ifdef G4HEPREPFILEDEBUG
	G4cout <<
    "G4HepRepFileSceneHandler::AddSolid(const G4Box& box) called for "
	<< box.GetName()
	<< G4endl;
	PrintThings();
#endif
	
	if (drawingTraj)
		return;
	
	if (drawingHit)
		InitHit();
	
	haveVisible = false;
	AddHepRepInstance("Prism", NULL);
	
	G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();
	
	// Get and check applicable vis attributes.
	fpVisAttribs = fpViewer->GetApplicableVisAttributes(fpVisAttribs);
	if (fpVisAttribs && (fpVisAttribs->IsVisible()==0) && messenger->getCullInvisibles())
		return;
	
	hepRepXMLWriter->addPrimitive();
	
	G4double dx = box.GetXHalfLength();
	G4double dy = box.GetYHalfLength();
	G4double dz = box.GetZHalfLength();
	
	G4Point3D vertex1(G4Point3D( dx, dy,-dz));
	G4Point3D vertex2(G4Point3D( dx,-dy,-dz));
	G4Point3D vertex3(G4Point3D(-dx,-dy,-dz));
	G4Point3D vertex4(G4Point3D(-dx, dy,-dz));
	G4Point3D vertex5(G4Point3D( dx, dy, dz));
	G4Point3D vertex6(G4Point3D( dx,-dy, dz));
	G4Point3D vertex7(G4Point3D(-dx,-dy, dz));
	G4Point3D vertex8(G4Point3D(-dx, dy, dz));
	
	vertex1 = (fObjectTransformation) * vertex1;
	vertex2 = (fObjectTransformation) * vertex2;
	vertex3 = (fObjectTransformation) * vertex3;
	vertex4 = (fObjectTransformation) * vertex4;
	vertex5 = (fObjectTransformation) * vertex5;
	vertex6 = (fObjectTransformation) * vertex6;
	vertex7 = (fObjectTransformation) * vertex7;
	vertex8 = (fObjectTransformation) * vertex8;
	
	hepRepXMLWriter->addPoint(vertex1.x(), vertex1.y(), vertex1.z());
	hepRepXMLWriter->addPoint(vertex2.x(), vertex2.y(), vertex2.z());
	hepRepXMLWriter->addPoint(vertex3.x(), vertex3.y(), vertex3.z());
	hepRepXMLWriter->addPoint(vertex4.x(), vertex4.y(), vertex4.z());
	hepRepXMLWriter->addPoint(vertex5.x(), vertex5.y(), vertex5.z());
	hepRepXMLWriter->addPoint(vertex6.x(), vertex6.y(), vertex6.z());
	hepRepXMLWriter->addPoint(vertex7.x(), vertex7.y(), vertex7.z());
	hepRepXMLWriter->addPoint(vertex8.x(), vertex8.y(), vertex8.z());
}


void G4HepRepFileSceneHandler::AddSolid(const G4Cons& cons) {
#ifdef G4HEPREPFILEDEBUG
	G4cout <<
    "G4HepRepFileSceneHandler::AddSolid(const G4Cons& cons) called for "
	<< cons.GetName()
	<< G4endl;
	PrintThings();
#endif
	
	// HepRApp does not correctly represent the end faces of cones at
	// non-standard angles, let the base class convert these solids to polygons.	
	G4RotationMatrix r = fObjectTransformation.getRotation();	
	G4bool linedUpWithAnAxis = (std::fabs(r.phiX())<=.001 ||  
								std::fabs(r.phiY())<=.001 || 
								std::fabs(r.phiZ())<=.001 ||
								std::fabs(r.phiX()-pi)<=.001 ||
								std::fabs(r.phiY()-pi)<=.001 ||
								std::fabs(r.phiZ()-pi)<=.001);	
	//G4cout << "Angle X:" << r.phiX() << ", Angle Y:" << r.phiY() << ", Angle Z:" << r.phiZ() << G4endl;
	//G4cout << "linedUpWithAnAxis:" << linedUpWithAnAxis << G4endl;
	
	// HepRep does not have a primitive for a cut cone,
	// so if this cone is cut, let the base class convert this
	// solid to polygons.
	G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();
	if (cons.GetDeltaPhiAngle() < twopi || !linedUpWithAnAxis || messenger->renderCylAsPolygons())
	{
		G4VSceneHandler::AddSolid(cons);  // Invoke default action.
	} else {
	
		if (drawingTraj)
			return;
		
		if (drawingHit)
			InitHit();

		haveVisible = false;
		AddHepRepInstance("Cylinder", NULL);
		
		// Get and check applicable vis attributes.
		fpVisAttribs = fpViewer->GetApplicableVisAttributes(fpVisAttribs);
		if (fpVisAttribs && (fpVisAttribs->IsVisible()==0) && messenger->getCullInvisibles())
			return;
		
		G4Point3D vertex1(G4Point3D( 0., 0., -cons.GetZHalfLength()));
		G4Point3D vertex2(G4Point3D( 0., 0.,  cons.GetZHalfLength()));
		
		vertex1 = (fObjectTransformation) * vertex1;
		vertex2 = (fObjectTransformation) * vertex2;
				
		// Outer cylinder.
		hepRepXMLWriter->addPrimitive();
		hepRepXMLWriter->addAttValue("Radius1",messenger->getScale() * cons.GetOuterRadiusMinusZ());
		hepRepXMLWriter->addAttValue("Radius2",messenger->getScale() * cons.GetOuterRadiusPlusZ());
		hepRepXMLWriter->addPoint(vertex1.x(), vertex1.y(), vertex1.z());
		hepRepXMLWriter->addPoint(vertex2.x(), vertex2.y(), vertex2.z());
		
		// Inner cylinder.
		hepRepXMLWriter->addPrimitive();
		hepRepXMLWriter->addAttValue("Radius1",messenger->getScale() * cons.GetInnerRadiusMinusZ());
		hepRepXMLWriter->addAttValue("Radius2",messenger->getScale() * cons.GetInnerRadiusPlusZ());
		hepRepXMLWriter->addPoint(vertex1.x(), vertex1.y(), vertex1.z());
		hepRepXMLWriter->addPoint(vertex2.x(), vertex2.y(), vertex2.z());
	}
}


void G4HepRepFileSceneHandler::AddSolid(const G4Tubs& tubs) {
#ifdef G4HEPREPFILEDEBUG
	G4cout <<
    "G4HepRepFileSceneHandler::AddSolid(const G4Tubs& tubs) called for "
	<< tubs.GetName()
	<< G4endl;
	PrintThings();
#endif
	
	// HepRApp does not correctly represent the end faces of cylinders at
	// non-standard angles, let the base class convert these solids to polygons.	
	G4RotationMatrix r = fObjectTransformation.getRotation();	
	G4bool linedUpWithAnAxis = (std::fabs(r.phiX())<=.001 ||  
								std::fabs(r.phiY())<=.001 || 
								std::fabs(r.phiZ())<=.001 ||
								std::fabs(r.phiX()-pi)<=.001 ||
								std::fabs(r.phiY()-pi)<=.001 ||
								std::fabs(r.phiZ()-pi)<=.001);	
	//G4cout << "Angle X:" << r.phiX() << ", Angle Y:" << r.phiY() << ", Angle Z:" << r.phiZ() << G4endl;
	//G4cout << "linedUpWithAnAxis:" << linedUpWithAnAxis << G4endl;
	
	// HepRep does not have a primitive for a cut cylinder,
	// so if this cylinder is cut, let the base class convert this
	// solid to polygons.
	G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();
	if (tubs.GetDeltaPhiAngle() < twopi || !linedUpWithAnAxis || messenger->renderCylAsPolygons())
	{
		G4VSceneHandler::AddSolid(tubs);  // Invoke default action.
	} else {
	
		if (drawingTraj)
			return;
		
		if (drawingHit)
			InitHit();

		haveVisible = false;
		AddHepRepInstance("Cylinder", NULL);
		
		// Get and check applicable vis attributes.
		fpVisAttribs = fpViewer->GetApplicableVisAttributes(fpVisAttribs);
		if (fpVisAttribs && (fpVisAttribs->IsVisible()==0) && messenger->getCullInvisibles())
			return;
		
		G4Point3D vertex1(G4Point3D( 0., 0., -tubs.GetZHalfLength()));
		G4Point3D vertex2(G4Point3D( 0., 0.,  tubs.GetZHalfLength()));
		
		vertex1 = (fObjectTransformation) * vertex1;
		vertex2 = (fObjectTransformation) * vertex2;
		
		// Outer cylinder.
		hepRepXMLWriter->addPrimitive();
		hepRepXMLWriter->addAttValue("Radius1", messenger->getScale() * tubs.GetOuterRadius());
		hepRepXMLWriter->addAttValue("Radius2", messenger->getScale() * tubs.GetOuterRadius());
		hepRepXMLWriter->addPoint(vertex1.x(), vertex1.y(), vertex1.z());
		hepRepXMLWriter->addPoint(vertex2.x(), vertex2.y(), vertex2.z());
		
		// Inner cylinder.
		if (tubs.GetInnerRadius() != 0.) {
			hepRepXMLWriter->addPrimitive();
			hepRepXMLWriter->addAttValue("Radius1", messenger->getScale() * tubs.GetInnerRadius());
			hepRepXMLWriter->addAttValue("Radius2", messenger->getScale() * tubs.GetInnerRadius());
			hepRepXMLWriter->addPoint(vertex1.x(), vertex1.y(), vertex1.z());
			hepRepXMLWriter->addPoint(vertex2.x(), vertex2.y(), vertex2.z());
		}
	}
}


void G4HepRepFileSceneHandler::AddSolid(const G4Trd& trd) {
#ifdef G4HEPREPFILEDEBUG
	G4cout <<
    "G4HepRepFileSceneHandler::AddSolid(const G4Trd& trd) called for "
	<< trd.GetName()
	<< G4endl;
	PrintThings();
#endif
	
	if (drawingTraj)
		return;
	
	if (drawingHit)
		InitHit();
	
	haveVisible = false;
	AddHepRepInstance("Prism", NULL);
	
	G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();
	
	// Get and check applicable vis attributes.
	fpVisAttribs = fpViewer->GetApplicableVisAttributes(fpVisAttribs);
	if (fpVisAttribs && (fpVisAttribs->IsVisible()==0) && messenger->getCullInvisibles())
		return;
	
	hepRepXMLWriter->addPrimitive();
	
	G4double dx1 = trd.GetXHalfLength1();
	G4double dy1 = trd.GetYHalfLength1();
	G4double dx2 = trd.GetXHalfLength2();
	G4double dy2 = trd.GetYHalfLength2();
	G4double dz = trd.GetZHalfLength();
	
	G4Point3D vertex1(G4Point3D( dx1, dy1,-dz));
	G4Point3D vertex2(G4Point3D( dx1,-dy1,-dz));
	G4Point3D vertex3(G4Point3D(-dx1,-dy1,-dz));
	G4Point3D vertex4(G4Point3D(-dx1, dy1,-dz));
	G4Point3D vertex5(G4Point3D( dx2, dy2, dz));
	G4Point3D vertex6(G4Point3D( dx2,-dy2, dz));
	G4Point3D vertex7(G4Point3D(-dx2,-dy2, dz));
	G4Point3D vertex8(G4Point3D(-dx2, dy2, dz));
	
	vertex1 = (fObjectTransformation) * vertex1;
	vertex2 = (fObjectTransformation) * vertex2;
	vertex3 = (fObjectTransformation) * vertex3;
	vertex4 = (fObjectTransformation) * vertex4;
	vertex5 = (fObjectTransformation) * vertex5;
	vertex6 = (fObjectTransformation) * vertex6;
	vertex7 = (fObjectTransformation) * vertex7;
	vertex8 = (fObjectTransformation) * vertex8;
	
	hepRepXMLWriter->addPoint(vertex1.x(), vertex1.y(), vertex1.z());
	hepRepXMLWriter->addPoint(vertex2.x(), vertex2.y(), vertex2.z());
	hepRepXMLWriter->addPoint(vertex3.x(), vertex3.y(), vertex3.z());
	hepRepXMLWriter->addPoint(vertex4.x(), vertex4.y(), vertex4.z());
	hepRepXMLWriter->addPoint(vertex5.x(), vertex5.y(), vertex5.z());
	hepRepXMLWriter->addPoint(vertex6.x(), vertex6.y(), vertex6.z());
	hepRepXMLWriter->addPoint(vertex7.x(), vertex7.y(), vertex7.z());
	hepRepXMLWriter->addPoint(vertex8.x(), vertex8.y(), vertex8.z());
}


void G4HepRepFileSceneHandler::AddSolid(const G4Trap& trap) {
#ifdef G4HEPREPFILEDEBUG
	G4cout <<
    "G4HepRepFileSceneHandler::AddSolid(const G4Trap& trap) called for "
	<< trap.GetName()
	<< G4endl;
	PrintThings();
#endif
	G4VSceneHandler::AddSolid(trap);  // Invoke default action.
}


void G4HepRepFileSceneHandler::AddSolid(const G4Sphere& sphere) {
#ifdef G4HEPREPFILEDEBUG
	G4cout <<
    "G4HepRepFileSceneHandler::AddSolid(const G4Sphere& sphere) called for "
	<< sphere.GetName()
	<< G4endl;
	PrintThings();
#endif
	G4VSceneHandler::AddSolid(sphere);  // Invoke default action.
}


void G4HepRepFileSceneHandler::AddSolid(const G4Para& para) {
#ifdef G4HEPREPFILEDEBUG
	G4cout <<
    "G4HepRepFileSceneHandler::AddSolid(const G4Para& para) called for "
	<< para.GetName()
	<< G4endl;
	PrintThings();
#endif
	G4VSceneHandler::AddSolid(para);  // Invoke default action.
}


void G4HepRepFileSceneHandler::AddSolid(const G4Torus& torus) {
#ifdef G4HEPREPFILEDEBUG
	G4cout <<
    "G4HepRepFileSceneHandler::AddSolid(const G4Torus& torus) called for "
	<< torus.GetName()
	<< G4endl;
	PrintThings();
#endif
	G4VSceneHandler::AddSolid(torus);  // Invoke default action.
}


void G4HepRepFileSceneHandler::AddSolid(const G4Polycone& polycone) {
#ifdef G4HEPREPFILEDEBUG
  G4cout <<
  "G4HepRepFileSceneHandler::AddSolid(const G4Polycone& polycone) called for "
  << polycone.GetName()
  << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddSolid(polycone);  // Invoke default action.
}


void G4HepRepFileSceneHandler::AddSolid(const G4Polyhedra& polyhedra) {
#ifdef G4HEPREPFILEDEBUG
  G4cout <<
  "G4HepRepFileSceneHandler::AddSolid(const G4Polyhedra& polyhedra) called for "
  << polyhedra.GetName()
  << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddSolid(polyhedra);  // Invoke default action.
}


void G4HepRepFileSceneHandler::AddSolid(const G4Orb& orb) {
#ifdef G4HEPREPFILEDEBUG
  G4cout <<
  "G4HepRepFileSceneHandler::AddSolid(const G4Orb& orb) called for "
  << orb.GetName()
  << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddSolid(orb);  // Invoke default action.
}


void G4HepRepFileSceneHandler::AddSolid(const G4Ellipsoid& ellipsoid) {
#ifdef G4HEPREPFILEDEBUG
  G4cout <<
  "G4HepRepFileSceneHandler::AddSolid(const G4Ellipsoid& ellipsoid) called for "
  << ellipsoid.GetName()
  << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddSolid(ellipsoid);  // Invoke default action.
}


void G4HepRepFileSceneHandler::AddSolid(const G4VSolid& solid) {
#ifdef G4HEPREPFILEDEBUG
	G4cout <<
    "G4HepRepFileSceneHandler::AddSolid(const G4Solid& solid) called for "
	<< solid.GetName()
	<< G4endl;
	PrintThings();
#endif
	G4VSceneHandler::AddSolid(solid);  // Invoke default action.
}


void G4HepRepFileSceneHandler::AddCompound (const G4VTrajectory& traj) {
#ifdef G4HEPREPFILEDEBUG
	G4cout << "G4HepRepFileSceneHandler::AddCompound(const G4VTrajectory&) " << G4endl;
#endif
	
	G4TrajectoriesModel* pTrModel =
		dynamic_cast<G4TrajectoriesModel*>(fpModel);
	if (!pTrModel) G4Exception
	  ("G4HepRepFileSceneHandler::AddCompound(const G4VTrajectory&)",
	   "vis-HepRep0001", FatalException, "Not a G4TrajectoriesModel.");
	
	// Pointers to hold trajectory attribute values and definitions.
	std::vector<G4AttValue>* rawTrajAttValues = traj.CreateAttValues();
	trajAttValues =
		new std::vector<G4AttValue>;
	trajAttDefs =
		new std::map<G4String,G4AttDef>;
	
	// Iterators to use with attribute values and definitions.
	std::vector<G4AttValue>::iterator iAttVal;
	std::map<G4String,G4AttDef>::const_iterator iAttDef;
	G4int i;
	
	// Get trajectory attributes and definitions in standard HepRep style
	// (uniform units, 3Vectors decomposed).
	if (rawTrajAttValues) {
		G4bool error = G4AttCheck(rawTrajAttValues,
								  traj.GetAttDefs()).Standard(trajAttValues,trajAttDefs);
		if (error) {
			G4cout << "G4HepRepFileSceneHandler::AddCompound(traj):"
			"\nERROR found during conversion to standard trajectory attributes."
			<< G4endl;
		}
#ifdef G4HEPREPFILEDEBUG 
		G4cout <<
			"G4HepRepFileSceneHandler::AddCompound(traj): standardised attributes:\n"
			<< G4AttCheck(trajAttValues,trajAttDefs) << G4endl;
#endif
		delete rawTrajAttValues;
	}
	
	// Open the HepRep output file if it is not already open.
	CheckFileOpen();
	
	// Add the Event Data Type if it hasn't already been added.
	if (strcmp("Event Data",hepRepXMLWriter->prevTypeName[0])!=0) {
		hepRepXMLWriter->addType("Event Data",0);
		hepRepXMLWriter->addInstance();
	}
	
	// Add the Trajectories Type.
	G4String previousName = hepRepXMLWriter->prevTypeName[1];
	hepRepXMLWriter->addType("Trajectories",1);
	
	// If this is the first trajectory of this event,
	// specify attribute values common to all trajectories.
	if (strcmp("Trajectories",previousName)!=0) {
		hepRepXMLWriter->addAttValue("Layer",100);
		
		// Take all Trajectory attDefs from first trajectory.
		// Would rather be able to get these attDefs without needing a reference from any
		// particular trajectory, but don't know how to do that.
		// Write out trajectory attribute definitions.
		if (trajAttValues && trajAttDefs) {
			for (iAttVal = trajAttValues->begin();
				 iAttVal != trajAttValues->end(); ++iAttVal) {
				iAttDef = trajAttDefs->find(iAttVal->GetName());
				if (iAttDef != trajAttDefs->end()) {
					// Protect against incorrect use of Category.  Anything value other than the
					// standard ones will be considered to be in the physics category.
					G4String category = iAttDef->second.GetCategory();
					if (strcmp(category,"Draw")!=0 &&
						strcmp(category,"Physics")!=0 &&
						strcmp(category,"Association")!=0 &&
						strcmp(category,"PickAction")!=0)
						category = "Physics";
					hepRepXMLWriter->addAttDef(iAttVal->GetName(), iAttDef->second.GetDesc(),
											   category, iAttDef->second.GetExtra());
				}
			} 
		}
		
		// Take all TrajectoryPoint attDefs from first point of first trajectory.
		// Would rather be able to get these attDefs without needing a reference from any
		// particular point, but don't know how to do that.
		if ((trajContext->GetDrawStepPts() || trajContext->GetDrawAuxPts())
			&& traj.GetPointEntries()>0) {
			G4VTrajectoryPoint* aTrajectoryPoint = traj.GetPoint(0);
			
			// Pointers to hold trajectory point attribute values and definitions.
			std::vector<G4AttValue>* rawPointAttValues = aTrajectoryPoint->CreateAttValues();
			std::vector<G4AttValue>* pointAttValues =
			new std::vector<G4AttValue>;
			std::map<G4String,G4AttDef>* pointAttDefs =
			new std::map<G4String,G4AttDef>;
			
			// Get first trajectory point's attributes and definitions in standard HepRep style
			// (uniform units, 3Vectors decomposed).
			if (rawPointAttValues) {
				G4bool error = G4AttCheck(rawPointAttValues,
										  aTrajectoryPoint->GetAttDefs()).Standard(pointAttValues,pointAttDefs);
				if (error) {
					G4cout << "G4HepRepFileSceneHandler::AddCompound(traj):"
					"\nERROR found during conversion to standard first point attributes." << G4endl;
				}
				
				// Write out point attribute definitions.
				if (pointAttValues && pointAttDefs) {
					for (iAttVal = pointAttValues->begin();
						 iAttVal != pointAttValues->end(); ++iAttVal) {
						iAttDef =
						pointAttDefs->find(iAttVal->GetName());
						if (iAttDef != pointAttDefs->end()) {
							// Protect against incorrect use of Category.  Anything value other than the
							// standard ones will be considered to be in the physics category.
							G4String category = iAttDef->second.GetCategory();
							if (strcmp(category,"Draw")!=0 &&
								strcmp(category,"Physics")!=0 &&
								strcmp(category,"Association")!=0 &&
								strcmp(category,"PickAction")!=0)
								category = "Physics";
							// Do not write out the Aux or Pos attribute.  Aux does not conform to the HepRep rule
							// that each object can have only one instance of a given AttValue.
							// Both of these attributes are redundant to actual position information of the point.
							if (strcmp(iAttVal->GetName(),"Aux-X")!=0 &&
								strcmp(iAttVal->GetName(),"Aux-Y")!=0 &&
								strcmp(iAttVal->GetName(),"Aux-Z")!=0 &&
								strcmp(iAttVal->GetName(),"Pos-X")!=0 &&
								strcmp(iAttVal->GetName(),"Pos-Y")!=0 &&
								strcmp(iAttVal->GetName(),"Pos-Z")!=0)
								hepRepXMLWriter->addAttDef(iAttVal->GetName(), iAttDef->second.GetDesc(),
														category, iAttDef->second.GetExtra());
						}
					}
				}
				delete rawPointAttValues;
			}
			
			// Clean up point attributes.
			if (pointAttValues)
				delete pointAttValues;
			if (pointAttDefs)
				delete pointAttDefs;
		}
	} // end of special treatment for when this is the first trajectory.
	
	// Now that we have written out all of the attributes that are based on the
	// trajectory's particulars, call base class to deconstruct trajectory into polyline and/or points
	// (or nothing if trajectory is to be filtered out).
	// If base class calls for drawing points, no points will actually be drawn there since we
	// instead need to do point drawing from here (in order to obtain the points attributes,
	// not available from AddPrimitive(...point).  Instead, such a call will just serve to set the
	// flag that tells us that point drawing was requested for this trajectory (depends on several
	// factors including trajContext and filtering).
	drawingTraj = true;
	doneInitTraj = false;
	G4VSceneHandler::AddCompound(traj);
	drawingTraj = false;
	
	// Draw step points.
	if (trajContext->GetDrawStepPts()) {
		if (!doneInitTraj)
			InitTrajectory();
		// Create Trajectory Points as a subType of Trajectories.
		// Note that we should create this heprep type even if there are no actual points.
		// This allows the user to tell that points don't exist (admittedly odd) rather
		// than that they were omitted by the drawing mode.
		previousName = hepRepXMLWriter->prevTypeName[2];
		hepRepXMLWriter->addType("Trajectory Step Points",2);

		float redness;
		float greenness;
		float blueness;			
		G4int markSize;
		G4bool visible;
		G4bool square;
		G4Colour colour = trajContext->GetStepPtsColour();
		redness = colour.GetRed();
		greenness = colour.GetGreen();
		blueness = colour.GetBlue();
		markSize = (G4int) trajContext->GetStepPtsSize();
		visible = (G4int) trajContext->GetStepPtsVisible();
		square = (trajContext->GetStepPtsType()==G4Polymarker::squares);

		// Avoiding drawing anything black on black.  
		if (redness==0. && greenness==0. && blueness==0.) {
			redness = 1.;
			greenness = 1.;
			blueness = 1.;
		}
		
		// Specify attributes common to all trajectory points.
		if (strcmp("Trajectory Step Points",previousName)!=0) {
			hepRepXMLWriter->addAttValue("DrawAs","Point");
			hepRepXMLWriter->addAttValue("MarkColor", redness, greenness, blueness);
			hepRepXMLWriter->addAttValue("MarkSize",markSize);
			hepRepXMLWriter->addAttValue("Layer",110);
			hepRepXMLWriter->addAttValue("Visibility",visible);
			if (square)
				hepRepXMLWriter->addAttValue("MarkName","square");
			else
				hepRepXMLWriter->addAttValue("MarkName","dot");
		}
		
		// Loop over all points on this trajectory.
		for (i = 0; i < traj.GetPointEntries(); i++) {
			G4VTrajectoryPoint* aTrajectoryPoint = traj.GetPoint(i);
			
			// Each point is a separate instance of the type Trajectory Points.
			hepRepXMLWriter->addInstance();
			
			// Pointers to hold trajectory point attribute values and definitions.
			std::vector<G4AttValue>* rawPointAttValues = aTrajectoryPoint->CreateAttValues();
			std::vector<G4AttValue>* pointAttValues =
				new std::vector<G4AttValue>;
			std::map<G4String,G4AttDef>* pointAttDefs =
				new std::map<G4String,G4AttDef>;
			
			// Get trajectory point attributes and definitions in standard HepRep style
			// (uniform units, 3Vectors decomposed).
			if (rawPointAttValues) {
				G4bool error = G4AttCheck(rawPointAttValues,
										  aTrajectoryPoint->GetAttDefs()).Standard(pointAttValues,pointAttDefs);
				if (error) {
					G4cout << "G4HepRepFileSceneHandler::AddCompound(traj):"
					"\nERROR found during conversion to standard point attributes." << G4endl;
				}
				
				// Write out point attribute values.
				if (pointAttValues) {
					for (iAttVal = pointAttValues->begin();
						 iAttVal != pointAttValues->end(); ++iAttVal)
						// Do not write out the Aux or Pos attribute.  Aux does not conform to the HepRep rule
						// that each object can have only one instance of a given AttValue.
						// Both of these attributes are redundant to actual position information of the point.
						if (strcmp(iAttVal->GetName(),"Aux-X")!=0 &&
							strcmp(iAttVal->GetName(),"Aux-Y")!=0 &&
							strcmp(iAttVal->GetName(),"Aux-Z")!=0 &&
							strcmp(iAttVal->GetName(),"Pos-X")!=0 &&
							strcmp(iAttVal->GetName(),"Pos-Y")!=0 &&
							strcmp(iAttVal->GetName(),"Pos-Z")!=0)
							hepRepXMLWriter->addAttValue(iAttVal->GetName(), iAttVal->GetValue());
				}
			}
			
			// Clean up point attributes.
                        delete pointAttDefs;
                        delete pointAttValues;
                        delete rawPointAttValues;

			// Each trajectory point is made of a single primitive, a point.
			hepRepXMLWriter->addPrimitive();
			G4Point3D vertex = aTrajectoryPoint->GetPosition();
			hepRepXMLWriter->addPoint(vertex.x(), vertex.y(), vertex.z());
		}
	}
	
	// Draw Auxiliary Points
	if (trajContext->GetDrawAuxPts()) {
		if (!doneInitTraj)
			InitTrajectory();
		// Create Trajectory Points as a subType of Trajectories.
		// Note that we should create this heprep type even if there are no actual points.
		// This allows the user to tell that points don't exist (admittedly odd) rather
		// than that they were omitted by the drawing mode.
		previousName = hepRepXMLWriter->prevTypeName[2];
		hepRepXMLWriter->addType("Trajectory Auxiliary Points",2);

		float redness;
		float greenness;
		float blueness;			
		G4int markSize;
		G4bool visible;
		G4bool square;
		G4Colour colour = trajContext->GetAuxPtsColour();
		redness = colour.GetRed();
		greenness = colour.GetGreen();
		blueness = colour.GetBlue();
		markSize = (G4int) trajContext->GetAuxPtsSize();
		visible = (G4int) trajContext->GetAuxPtsVisible();
		square = (trajContext->GetAuxPtsType()==G4Polymarker::squares);

		// Avoiding drawing anything black on black.  
		if (redness==0. && greenness==0. && blueness==0.) {
			redness = 1.;
			greenness = 1.;
			blueness = 1.;
		}
		
		// Specify attributes common to all trajectory points.
		if (strcmp("Trajectory Auxiliary Points",previousName)!=0) {
			hepRepXMLWriter->addAttValue("DrawAs","Point");
			hepRepXMLWriter->addAttValue("MarkColor", redness, greenness, blueness);
			hepRepXMLWriter->addAttValue("MarkSize",markSize);
			hepRepXMLWriter->addAttValue("Layer",110);
			hepRepXMLWriter->addAttValue("Visibility",visible);
			if (square)
				hepRepXMLWriter->addAttValue("MarkName","Square");
			else
				hepRepXMLWriter->addAttValue("MarkName","Dot");
		}
		
		// Loop over all points on this trajectory.
		for (i = 0; i < traj.GetPointEntries(); i++) {
			G4VTrajectoryPoint* aTrajectoryPoint = traj.GetPoint(i);
			
			// Each point is a separate instance of the type Trajectory Points.
			hepRepXMLWriter->addInstance();
			
			// Pointers to hold trajectory point attribute values and definitions.
			std::vector<G4AttValue>* rawPointAttValues = aTrajectoryPoint->CreateAttValues();
			std::vector<G4AttValue>* pointAttValues =
				new std::vector<G4AttValue>;
			std::map<G4String,G4AttDef>* pointAttDefs =
				new std::map<G4String,G4AttDef>;
			
			// Get trajectory point attributes and definitions in standard HepRep style
			// (uniform units, 3Vectors decomposed).
			if (rawPointAttValues) {
				G4bool error = G4AttCheck(rawPointAttValues,
										  aTrajectoryPoint->GetAttDefs()).Standard(pointAttValues,pointAttDefs);
				if (error) {
					G4cout << "G4HepRepFileSceneHandler::AddCompound(traj):"
					"\nERROR found during conversion to standard point attributes." << G4endl;
				}
				
				// Write out point attribute values.
				if (pointAttValues) {
					for (iAttVal = pointAttValues->begin();
						 iAttVal != pointAttValues->end(); ++iAttVal)
						// Do not write out the Aux or Pos attribute.  Aux does not conform to the HepRep rule
						// that each object can have only one instance of a given AttValue.
						// Both of these attributes are redundant to actual position information of the point.
						if (strcmp(iAttVal->GetName(),"Aux-X")!=0 &&
							strcmp(iAttVal->GetName(),"Aux-Y")!=0 &&
							strcmp(iAttVal->GetName(),"Aux-Z")!=0 &&
							strcmp(iAttVal->GetName(),"Pos-X")!=0 &&
							strcmp(iAttVal->GetName(),"Pos-Y")!=0 &&
							strcmp(iAttVal->GetName(),"Pos-Z")!=0)
							hepRepXMLWriter->addAttValue(iAttVal->GetName(), iAttVal->GetValue());
				}
			}
			
			// Clean up point attributes.
			delete pointAttDefs;
                        delete pointAttValues;
                        delete rawPointAttValues;

			// Each trajectory point is made of a single primitive, a point.
			G4Point3D vertex = aTrajectoryPoint->GetPosition();
			
			// Loop over auxiliary points associated with this Trajectory Point.
			const std::vector<G4ThreeVector>* auxiliaries = aTrajectoryPoint->GetAuxiliaryPoints();
			if (0 != auxiliaries) {
				for (size_t iAux=0; iAux<auxiliaries->size(); ++iAux) {
					const G4ThreeVector auxPos((*auxiliaries)[iAux]);
					hepRepXMLWriter->addPrimitive();
					hepRepXMLWriter->addPoint(auxPos.x(), auxPos.y(), auxPos.z());
				}
			}
		}
	}
}


void G4HepRepFileSceneHandler::AddCompound (const G4VHit& hit) {
#ifdef G4HEPREPFILEDEBUG
	G4cout << "G4HepRepFileSceneHandler::AddCompound(G4VHit&) " << G4endl;
#endif
	
	// Pointers to hold hit attribute values and definitions.
	std::vector<G4AttValue>* rawHitAttValues = hit.CreateAttValues();
	hitAttValues =
		new std::vector<G4AttValue>;
	hitAttDefs =
		new std::map<G4String,G4AttDef>;
	
	// Iterators to use with attribute values and definitions.
	std::vector<G4AttValue>::iterator iAttVal;
	std::map<G4String,G4AttDef>::const_iterator iAttDef;
	
	// Get hit attributes and definitions in standard HepRep style
	// (uniform units, 3Vectors decomposed).
	if (rawHitAttValues) {
		G4bool error = G4AttCheck(rawHitAttValues,
								  hit.GetAttDefs()).Standard(hitAttValues,hitAttDefs);
		if (error) {
			G4cout << "G4HepRepFileSceneHandler::AddCompound(hit):"
			"\nERROR found during conversion to standard hit attributes."
			<< G4endl;
		}
#ifdef G4HEPREPFILEDEBUG 
		G4cout <<
			"G4HepRepFileSceneHandler::AddCompound(hit): standardised attributes:\n"
			<< G4AttCheck(hitAttValues,hitAttDefs) << G4endl;
#endif
		delete rawHitAttValues;
	}
	
	// Open the HepRep output file if it is not already open.
	CheckFileOpen();
	
	// Add the Event Data Type if it hasn't already been added.
	if (strcmp("Event Data",hepRepXMLWriter->prevTypeName[0])!=0) {
		hepRepXMLWriter->addType("Event Data",0);
		hepRepXMLWriter->addInstance();
	}
	
	// Find out the current HitType.
	G4String hitType = "Hits";
	if (hitAttValues) {
		G4bool found = false;
		for (iAttVal = hitAttValues->begin();
			 iAttVal != hitAttValues->end() && !found; ++iAttVal) {
			if (strcmp(iAttVal->GetName(),"HitType")==0) {
				hitType = iAttVal->GetValue();
				found = true;
			}
		}
	}
	
	// Add the Hits Type.
	G4String previousName = hepRepXMLWriter->prevTypeName[1];
	hepRepXMLWriter->addType(hitType,1);
	
	// If this is the first hit of this event,
	// specify attribute values common to all hits.
	if (strcmp(hitType,previousName)!=0) {
		hepRepXMLWriter->addAttValue("Layer",130);
		
		// Take all Hit attDefs from first hit.
		// Would rather be able to get these attDefs without needing a reference from any
		// particular hit, but don't know how to do that.
		// Write out hit attribute definitions.
		if (hitAttValues && hitAttDefs) {
			for (iAttVal = hitAttValues->begin();
				 iAttVal != hitAttValues->end(); ++iAttVal) {
				iAttDef = hitAttDefs->find(iAttVal->GetName());
				if (iAttDef != hitAttDefs->end()) {
					// Protect against incorrect use of Category.  Anything value other than the
					// standard ones will be considered to be in the physics category.
					G4String category = iAttDef->second.GetCategory();
					if (strcmp(category,"Draw")!=0 &&
						strcmp(category,"Physics")!=0 &&
						strcmp(category,"Association")!=0 &&
						strcmp(category,"PickAction")!=0)
						category = "Physics";
					hepRepXMLWriter->addAttDef(iAttVal->GetName(), iAttDef->second.GetDesc(),
											   category, iAttDef->second.GetExtra());
				}
			}
		}
	} // end of special treatment for when this is the first hit.
	
	// Now that we have written out all of the attributes that are based on the
	// hit's particulars, call base class to deconstruct hit into a primitives.
	drawingHit = true;
	doneInitHit = false;
	G4VSceneHandler::AddCompound(hit);  // Invoke default action.
	drawingHit = false;
}


void G4HepRepFileSceneHandler::InitTrajectory() {
	if (!doneInitTraj) {
		// For every trajectory, add an instance of Type Trajectory.
		hepRepXMLWriter->addInstance();
		
		// Write out the trajectory's attribute values.
		if (trajAttValues) {
			std::vector<G4AttValue>::iterator iAttVal;
			for (iAttVal = trajAttValues->begin();
				 iAttVal != trajAttValues->end(); ++iAttVal)
				hepRepXMLWriter->addAttValue(iAttVal->GetName(), iAttVal->GetValue());
			delete trajAttValues; 
		}
		
		// Clean up trajectory attributes.
		if (trajAttDefs)
			delete trajAttDefs;
		
		doneInitTraj = true;
	}
}


void G4HepRepFileSceneHandler::InitHit() {
	if (!doneInitHit) {
		// For every hit, add an instance of Type Hit.
		hepRepXMLWriter->addInstance();
		
		// Write out the hit's attribute values.
		if (hitAttValues) {
			std::vector<G4AttValue>::iterator iAttVal;
			for (iAttVal = hitAttValues->begin();
				 iAttVal != hitAttValues->end(); ++iAttVal)
				hepRepXMLWriter->addAttValue(iAttVal->GetName(), iAttVal->GetValue());
			delete hitAttValues;
		}
		
		// Clean up hit attributes.
		if (hitAttDefs)
			delete hitAttDefs;
		
		doneInitHit = true;
	}
}


void G4HepRepFileSceneHandler::AddPrimitive(const G4Polyline& polyline) {
#ifdef G4HEPREPFILEDEBUG
	G4cout <<
    "G4HepRepFileSceneHandler::AddPrimitive(const G4Polyline& polyline) called:"
    "\n  polyline: " << polyline
	<< G4endl;
	PrintThings();
#endif
	
	G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();
	
	if (fpVisAttribs && (fpVisAttribs->IsVisible()==0) && messenger->getCullInvisibles())
		return;
	
	if (inPrimitives2D) {
		if (!warnedAbout2DMarkers) {
			G4cout << "HepRepFile does not currently support 2D lines." << G4endl;
			warnedAbout2DMarkers = true;
		}
		return;
    }
	
	if (drawingTraj)
		InitTrajectory();
		
	if (drawingHit)
		InitHit();
	
	haveVisible = true;
	AddHepRepInstance("Line", polyline);
	
	hepRepXMLWriter->addPrimitive();
	
	for (size_t i=0; i < polyline.size(); i++) {
		G4Point3D vertex = (fObjectTransformation) * polyline[i];
		hepRepXMLWriter->addPoint(vertex.x(), vertex.y(), vertex.z());
	}
}



void G4HepRepFileSceneHandler::AddPrimitive (const G4Polymarker& line) {
#ifdef G4HEPREPFILEDEBUG
	G4cout <<
    "G4HepRepFileSceneHandler::AddPrimitive(const G4Polymarker& line) called"
	<< G4endl;
	PrintThings();
#endif
	
	G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();
	
	if (fpVisAttribs && (fpVisAttribs->IsVisible()==0) && messenger->getCullInvisibles())
		return;
	
	if (inPrimitives2D) {
		if (!warnedAbout2DMarkers) {
			G4cout << "HepRepFile does not currently support 2D lines." << G4endl;
			warnedAbout2DMarkers = true;
		}
		return;
    }
	
	MarkerSizeType sizeType;
	G4double size = GetMarkerSize (line, sizeType);
	if (sizeType==world)
		size = 4.;
	
	if (drawingTraj)
		return;
	
	if (drawingHit)
		InitHit();
	
	haveVisible = true;
	AddHepRepInstance("Point", line);
	
	hepRepXMLWriter->addAttValue("MarkName", "Dot");
	hepRepXMLWriter->addAttValue("MarkSize", (G4int) size);
	
	hepRepXMLWriter->addPrimitive();
	
	for (size_t i=0; i < line.size(); i++) {
		G4Point3D vertex = (fObjectTransformation) * line[i];
		hepRepXMLWriter->addPoint(vertex.x(), vertex.y(), vertex.z());
	}
}


void G4HepRepFileSceneHandler::AddPrimitive(const G4Text& text) {
#ifdef G4HEPREPFILEDEBUG
	G4cout <<
    "G4HepRepFileSceneHandler::AddPrimitive(const G4Text& text) called:"
    "\n  text: " << text.GetText()
	<< G4endl;
	PrintThings();
#endif
	
	if (!inPrimitives2D) {
		if (!warnedAbout3DText) {
			G4cout << "HepRepFile does not currently support 3D text." << G4endl;
			G4cout << "HepRep browsers can directly display text attributes on request." << G4endl;
			G4cout << "See Application Developers Guide for how to attach attributes to viewable objects." << G4endl;
			warnedAbout3DText = true;
		}
		return;
    }
	
	MarkerSizeType sizeType;
	G4double size = GetMarkerSize (text, sizeType);
	if (sizeType==world)
		size = 12.;
	
	haveVisible = true;
	AddHepRepInstance("Text", text);
	
	hepRepXMLWriter->addAttValue("VAlign", "Top");
	hepRepXMLWriter->addAttValue("HAlign", "Left");
	hepRepXMLWriter->addAttValue("FontName", "Arial");
	hepRepXMLWriter->addAttValue("FontStyle", "Plain");
	hepRepXMLWriter->addAttValue("FontSize", (G4int) size);
	hepRepXMLWriter->addAttValue("FontHasBanner", "TRUE");
	hepRepXMLWriter->addAttValue("FontBannerColor", "0,0,0");
	
	const G4Colour& colour = GetTextColour(text);
	float redness = colour.GetRed();
	float greenness = colour.GetGreen();
	float blueness = colour.GetBlue();
    
	// Avoiding drawing anything black on black.  
	if (redness==0. && greenness==0. && blueness==0.) {
		redness = 1.;
		greenness = 1.;
		blueness = 1.;
	}
	hepRepXMLWriter->addAttValue("FontColor",redness,greenness,blueness);
	
	hepRepXMLWriter->addPrimitive();
	
	hepRepXMLWriter->addAttValue("Text", text.GetText());
	hepRepXMLWriter->addAttValue("VPos", .99-text.GetYOffset());
	hepRepXMLWriter->addAttValue("HPos", text.GetXOffset());
}


void G4HepRepFileSceneHandler::AddPrimitive(const G4Circle& circle) {
#ifdef G4HEPREPFILEDEBUG
	G4cout <<
    "G4HepRepFileSceneHandler::AddPrimitive(const G4Circle& circle) called:"
    "\n  radius: " << circle.GetWorldRadius()
	<< G4endl;
	PrintThings();
#endif
	
	G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();
	
	if (fpVisAttribs && (fpVisAttribs->IsVisible()==0) && messenger->getCullInvisibles())
		return;
	
	if (inPrimitives2D) {
		if (!warnedAbout2DMarkers) {
			G4cout << "HepRepFile does not currently support 2D circles." << G4endl;
			warnedAbout2DMarkers = true;
		}
		return;
    }
	
	MarkerSizeType sizeType;
	G4double size = GetMarkerSize (circle, sizeType);
	if (sizeType==world)
		size = 4.;
	
	if (drawingTraj)
		return;
	
	if (drawingHit)
		InitHit();
	
	haveVisible = true;
	AddHepRepInstance("Point", circle);
	
	hepRepXMLWriter->addAttValue("MarkName", "Dot");
	hepRepXMLWriter->addAttValue("MarkSize", (G4int) size);
	
	hepRepXMLWriter->addPrimitive();
	
	G4Point3D center = (fObjectTransformation) * circle.GetPosition();
	hepRepXMLWriter->addPoint(center.x(), center.y(), center.z());
}


void G4HepRepFileSceneHandler::AddPrimitive(const G4Square& square) {
#ifdef G4HEPREPFILEDEBUG
	G4cout <<
    "G4HepRepFileSceneHandler::AddPrimitive(const G4Square& square) called:"
    "\n  side: " << square.GetWorldRadius()
	<< G4endl;
	PrintThings();
#endif
	
	G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();
	
	if (fpVisAttribs && (fpVisAttribs->IsVisible()==0) && messenger->getCullInvisibles())
		return;
	
	if (inPrimitives2D) {
		if (!warnedAbout2DMarkers) {
			G4cout << "HepRepFile does not currently support 2D squares." << G4endl;
			warnedAbout2DMarkers = true;
		}
		return;
    }
	
	MarkerSizeType sizeType;
	G4double size = GetMarkerSize (square, sizeType);
	if (sizeType==world)
		size = 4.;
	
	if (drawingTraj)
		return;
	
	if (drawingHit)
		InitHit();
	
	haveVisible = true;
	AddHepRepInstance("Point", square);
	
	hepRepXMLWriter->addAttValue("MarkName", "Square");
	hepRepXMLWriter->addAttValue("MarkSize", (G4int) size);
	
	hepRepXMLWriter->addPrimitive();
	
	G4Point3D center = (fObjectTransformation) * square.GetPosition();
	hepRepXMLWriter->addPoint(center.x(), center.y(), center.z());
}


void G4HepRepFileSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) {
#ifdef G4HEPREPFILEDEBUG
	G4cout <<
    "G4HepRepFileSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) called."
	<< G4endl;
	PrintThings();
#endif
	
	G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();
	
	if (fpVisAttribs && (fpVisAttribs->IsVisible()==0) && messenger->getCullInvisibles())
		return;
	
	if(polyhedron.GetNoFacets()==0)return;
	
	if (drawingTraj)
		return;
	
	if (drawingHit)
		InitHit();
	
	haveVisible = true;
	AddHepRepInstance("Polygon", polyhedron);
	
	G4Normal3D surfaceNormal;
	G4Point3D vertex;
	
	G4bool notLastFace;
	do {
		hepRepXMLWriter->addPrimitive();
		notLastFace = polyhedron.GetNextNormal (surfaceNormal);
		
		G4int edgeFlag = 1;
		G4bool notLastEdge;
		do {
			notLastEdge = polyhedron.GetNextVertex (vertex, edgeFlag);
			vertex = (fObjectTransformation) * vertex;
			hepRepXMLWriter->addPoint(vertex.x(), vertex.y(), vertex.z());
		} while (notLastEdge);
	} while (notLastFace);
}


G4HepRepFileXMLWriter *G4HepRepFileSceneHandler::GetHepRepXMLWriter() {
    return hepRepXMLWriter;
}


void G4HepRepFileSceneHandler::AddHepRepInstance(const char* primName,
                                                 const G4Visible visible) {
#ifdef G4HEPREPFILEDEBUG
	G4cout <<
    "G4HepRepFileSceneHandler::AddHepRepInstance called."
	<< G4endl;
#endif
	
	// Open the HepRep output file if it is not already open.
	CheckFileOpen();
	
	G4VPhysicalVolume* pCurrentPV = 0;
	G4LogicalVolume* pCurrentLV = 0;
	G4int currentDepth = 0;
	G4PhysicalVolumeModel* pPVModel = dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
	if (pPVModel) {
		pCurrentPV = pPVModel->GetCurrentPV();
		pCurrentLV = pPVModel->GetCurrentLV();
		currentDepth = pPVModel->GetCurrentDepth();
	}
	
#ifdef G4HEPREPFILEDEBUG
	G4cout <<
		"pCurrentPV:" << pCurrentPV << ", readyForTransients:" << fReadyForTransients
		<< G4endl;
#endif
	
	if (drawingTraj || drawingHit) {
		// In this case, HepRep type, layer and instance were already created
		// in the AddCompound method.
	}
	else if (fReadyForTransients) {
		if (strcmp("Event Data",hepRepXMLWriter->prevTypeName[0])!=0) {
			hepRepXMLWriter->addType("Event Data",0);
			hepRepXMLWriter->addInstance();
		}
		
		// Applications have the option of either calling AddSolid(G4VTrajectory&) and
		// AddSolid(G4VHits&), or of just decomposing these into simpler primitives.
		// In the former case, drawing will be handled above and will include setting of
		// physics attributes.
		// In the latter case, which is an older style of working, we end up drawing the
		// trajectories and hits here, where we have no access to physics attributes. 
		// We receive primitives here.  We can figure out that these are transients, but we
		// have to guess exactly what these transients represent.
		// We assume the primitives are being used as in G4VTrajectory, hence we assume:
		// Lines are Trajectories
		// Squares that come after we've seen trajectories are Auxiliary Points
		// Circles that come after we've seen trajectories are Step Points
		// Other primitives are Hits
		
		int layer;
		
		if (strcmp("Text",primName)==0) {
			hepRepXMLWriter->addType("EventID",1);
		} else {
			if (strcmp("Line",primName)==0) {
				hepRepXMLWriter->addType("TransientPolylines",1);
				layer = 100;
			} else {
				if (strcmp(hepRepXMLWriter->prevTypeName[1],"TransientPolylines")==0 &&
					strcmp("Square",primName)==0)
				{
					hepRepXMLWriter->addType("AuxiliaryPoints",2);
					layer = 110;
				} else {
					if (strcmp(hepRepXMLWriter->prevTypeName[1],"TransientPolylines")==0 &&
						strcmp("Circle",primName)==0)
					{
						hepRepXMLWriter->addType("StepPoints",2);
						layer = 120;
					} else {
						hepRepXMLWriter->addType("Hits",1);
						layer = 130;
					}
				}
			}
			hepRepXMLWriter->addAttValue("Layer",layer);
		}
		
		hepRepXMLWriter->addInstance();
		
		// Handle Type declaration for Axes, Ruler, etc.
	} else if (pCurrentPV==0) {
		if (strcmp("AxesEtc",hepRepXMLWriter->prevTypeName[0])!=0) {
			hepRepXMLWriter->addType("AxesEtc",0);
			hepRepXMLWriter->addInstance();
		}
		
		int layer;
		
		if (strcmp("Text",primName)==0) {
			hepRepXMLWriter->addType("Text",1);
		} else {
			if (strcmp("Line",primName)==0) {
				hepRepXMLWriter->addType("Polylines",1);
				layer = 100;
			} else {
				hepRepXMLWriter->addType("Points",1);
				layer = 130;
			}
			hepRepXMLWriter->addAttValue("Layer",layer);
		}
		
		hepRepXMLWriter->addInstance();
		
		// Handle Type declaration for Detector Geometry,
		// replacing G4's top geometry level name "worldPhysical" with the
		// name "Detector Geometry".
	} else {
		//G4cout << "CurrentDepth" << currentDepth << G4endl;
		//G4cout << "currentName" << pCurrentPV->GetName() << G4endl;
		if (strcmp("Detector Geometry",hepRepXMLWriter->prevTypeName[0])!=0) {
			//G4cout << "Adding Det Geom type" << G4endl;
			hepRepXMLWriter->addType("Detector Geometry",0);
			hepRepXMLWriter->addInstance();
		}
		
		// Re-insert any layers of the hierarchy that were removed by G4's culling process.
		// Don't bother checking if same type name as last instance.
		if(strcmp(hepRepXMLWriter->prevTypeName[currentDepth+1],pCurrentPV->GetName())!=0) {
			//G4cout << "Looking for mother of:" << pCurrentLV->GetName() << G4endl;
			typedef G4PhysicalVolumeModel::G4PhysicalVolumeNodeID PVNodeID;
			typedef std::vector<PVNodeID> PVPath;
			const PVPath& drawnPVPath = pPVModel->GetDrawnPVPath();
			PVPath::const_reverse_iterator ri = ++drawnPVPath.rbegin();
			G4int drawnMotherDepth;
			if (ri != drawnPVPath.rend()) {
				// This volume has a mother.
				drawnMotherDepth = ri->GetNonCulledDepth();
				//G4cout << "drawnMotherDepth" << drawnMotherDepth << G4endl;
			} else {
				// This volume has no mother.  Must be a top level volume.
				drawnMotherDepth = -1;
				//G4cout << "Mother must be very top" << G4endl;
			}
			
			while (drawnMotherDepth < (currentDepth-1)) {
				G4String culledParentName = "Culled parent of " + pCurrentPV->GetName();
				//G4cout << "Inserting culled layer " << culledParentName << " at depth:" << drawnMotherDepth+2 << G4endl;
				hepRepXMLWriter->addType(culledParentName, drawnMotherDepth+2);
				hepRepXMLWriter->addInstance();
				drawnMotherDepth ++;
			}
		}
		
		// Add the HepRepType for the current volume.
		hepRepXMLWriter->addType(pCurrentPV->GetName(),currentDepth+1);
		hepRepXMLWriter->addInstance();
		
		G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();
		
		if (fpVisAttribs && (fpVisAttribs->IsVisible()==0) && messenger->getCullInvisibles())
			return;
		
		// Additional attributes.
		hepRepXMLWriter->addAttValue("Layer",hepRepXMLWriter->typeDepth);
		hepRepXMLWriter->addAttValue("LVol", pCurrentLV->GetName());
		G4Region* region = pCurrentLV->GetRegion();
		G4String regionName = region? region->GetName(): G4String("No region");
		hepRepXMLWriter->addAttValue("Region", regionName);
		hepRepXMLWriter->addAttValue("RootRegion", pCurrentLV->IsRootRegion());
		hepRepXMLWriter->addAttValue("Solid", pCurrentLV->GetSolid()->GetName());
		hepRepXMLWriter->addAttValue("EType", pCurrentLV->GetSolid()->GetEntityType());
		G4Material * material = pPVModel->GetCurrentMaterial();
		G4String matName = material? material->GetName(): G4String("No material");
		hepRepXMLWriter->addAttValue("Material", matName);
		G4double matDensity = material? material->GetDensity(): 0.;
		hepRepXMLWriter->addAttValue("Density", matDensity*m3/kg);
		G4State matState = material? material->GetState(): kStateUndefined;
		hepRepXMLWriter->addAttValue("State", matState);
		G4double matRadlen = material? material->GetRadlen(): 0.;
		hepRepXMLWriter->addAttValue("Radlen", matRadlen/m);
	}
	
	hepRepXMLWriter->addAttValue("DrawAs",primName);
	
	// Handle color and visibility attributes.
	float redness;
	float greenness;
	float blueness;
	G4bool isVisible; 
	
	if (fpVisAttribs || haveVisible) {
		G4Colour colour;
		
		if (fpVisAttribs) {
			colour = fpVisAttribs->GetColour();
			isVisible = fpVisAttribs->IsVisible();
		} else {
			colour = visible.GetVisAttributes()->GetColour();
			isVisible = fpViewer->
			GetApplicableVisAttributes(visible.GetVisAttributes())->IsVisible();
		}

		redness = colour.GetRed();
		greenness = colour.GetGreen();
		blueness = colour.GetBlue();
		
		// Avoiding drawing anything black on black.  
		if (redness==0. && greenness==0. && blueness==0.) {
			redness = 1.;
			greenness = 1.;
			blueness = 1.;
		}
	} else {
#ifdef G4HEPREPFILEDEBUG
		G4cout <<
		"G4HepRepFileSceneHandler::AddHepRepInstance using default colour."
		<< G4endl;
#endif
		redness = 1.;
		greenness = 1.;
		blueness = 1.;
		isVisible = true;
	}
	
	if (strcmp(primName,"Point")==0)
		hepRepXMLWriter->addAttValue("MarkColor",redness,greenness,blueness);
	else
		hepRepXMLWriter->addAttValue("LineColor",redness,greenness,blueness);
	
	hepRepXMLWriter->addAttValue("Visibility",isVisible);
}


void G4HepRepFileSceneHandler::CheckFileOpen() {
#ifdef G4HEPREPFILEDEBUG
	G4cout <<
    "G4HepRepFileSceneHandler::CheckFileOpen called."
	<< G4endl;
#endif
	
	if (!hepRepXMLWriter->isOpen) {
		G4String newFileSpec;

		G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();
		
		if (messenger->getOverwrite()) {
			newFileSpec = messenger->getFileDir()+messenger->getFileName()+".heprep";
		} else {
			newFileSpec = messenger->getFileDir()+messenger->getFileName()+G4UIcommand::ConvertToString(fileCounter)+".heprep";
		}
		
		G4cout << "HepRepFile writing to " << newFileSpec << G4endl;
		
		hepRepXMLWriter->open(newFileSpec);
		
		if (!messenger->getOverwrite())
			fileCounter++;
		
		hepRepXMLWriter->addAttDef("Generator", "HepRep Data Generator", "Physics","");
		G4String versionString = G4Version;
		versionString = versionString.substr(1,versionString.size()-2);
		versionString = " Geant4 version " + versionString + "   " + G4Date;
		hepRepXMLWriter->addAttValue("Generator", versionString);
		
		hepRepXMLWriter->addAttDef("LVol", "Logical Volume", "Physics","");
		hepRepXMLWriter->addAttDef("Region", "Cuts Region", "Physics","");
		hepRepXMLWriter->addAttDef("RootRegion", "Root Region", "Physics","");
		hepRepXMLWriter->addAttDef("Solid", "Solid Name", "Physics","");
		hepRepXMLWriter->addAttDef("EType", "Entity Type", "Physics","");
		hepRepXMLWriter->addAttDef("Material", "Material Name", "Physics","");
		hepRepXMLWriter->addAttDef("Density", "Material Density", "Physics","kg/m3");
		hepRepXMLWriter->addAttDef("State", "Material State", "Physics","");
		hepRepXMLWriter->addAttDef("Radlen", "Material Radiation Length", "Physics","m");
	}
}


void G4HepRepFileSceneHandler::ClearTransientStore() {
	// This is typically called after an update and before drawing hits
	// of the next event.  To simulate the clearing of "transients"
	// (hits, etc.) the detector is redrawn...
	if (fpViewer) {
		fpViewer -> SetView();
		fpViewer -> ClearView();
		fpViewer -> DrawView();
	}
}
