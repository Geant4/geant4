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
//

/**
 * @author Mark Donszelmann
 * @version $Id: G4HepRepSceneHandler.cc,v 1.8 2002-11-14 05:08:09 duns Exp $
 */

#include "g4std/vector"
#include "g4std/strstream"
#include "g4std/fstream"

//HepRep
#include "HEPREP/HepRep.h"
#include "HEPREP/HepRepAttribute.h"
#include "HEPREP/HepRepFactory.h"
#include "HEPREP/HepRepInstance.h"
#include "HEPREP/HepRepInstanceTree.h"
#include "HEPREP/HepRepPoint.h"
#include "HEPREP/HepRepTreeID.h"
#include "HEPREP/HepRepType.h"
#include "HEPREP/HepRepTypeTree.h"

//G4
#include "G4Types.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4Polyhedron.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Text.hh"
#include "G4NURBS.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4VSolid.hh"
#include "G4Scene.hh"

//This
#include "G4HepRep.hh"
#include "G4HepRepSceneHandler.hh"
#include "XMLHepRepStreamerFactory.h"


using namespace HEPREP;
using namespace std;

G4int G4HepRepSceneHandler::sceneCount = 0;
G4int G4HepRepSceneHandler::sceneIdCount = 0;

G4HepRepSceneHandler::G4HepRepSceneHandler (G4VGraphicsSystem& system, const G4String& name)
        : G4VSceneHandler (system, sceneIdCount++, name),
          currentDepth  (0),
          currentPV     (0),
          currentLV     (0) {

#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::G4HepRepSceneHandler: "
           << system << " "
           << GetScene()->GetName() << G4endl;
#endif
    heprepFactory = new XMLHepRepStreamerFactory();
    writer = NULL;
}

G4HepRepSceneHandler::~G4HepRepSceneHandler () {
    close();
    delete heprepFactory;
    heprepFactory = NULL;
}

HepRepFactory* G4HepRepSceneHandler::GetHepRepFactory() {
    if (writer == NULL) open();
    return heprepFactory;
}

void G4HepRepSceneHandler::open() {
    if (writer != NULL) return;

    char fname [256];
    ostrstream ost(fname, 256);
    ost << GetScene()->GetName() << ends;

#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::open(" << fname << ") " << G4endl;
#endif
    if (strcmp(fname, "stdout") == 0) {
        out = NULL;
        writer = heprepFactory->createHepRepWriter(&G4cout);
    } else if (strcmp(fname, "stderr") == 0) {
        out = NULL;
        writer = heprepFactory->createHepRepWriter(&G4cerr);
    } else {
        out = new ofstream(fname);
        writer = heprepFactory->createHepRepWriter(out);
    }

	heprep = heprepFactory->createHepRep();
	heprep->addLayer("Geometry");
	heprep->addLayer("Event");
	heprep->addLayer("CalHit");
	heprep->addLayer("Track");
	heprep->addLayer("Hit");

    HepRepTreeID* treeID = heprepFactory->createHepRepTreeID("G4Types", "1.0");
    HepRepTypeTree* typeTree = heprepFactory->createHepRepTypeTree(treeID);
    heprep->addTypeTree(typeTree);

    geometryType = heprepFactory->createHepRepType(NULL, "Detector Geometry");
    geometryType->addAttDef("LVol", "Logical Volume", "Physics","");
    geometryType->addAttDef("Solid", "Solid Name", "Physics","");
    geometryType->addAttDef("EType", "Entity Type", "Physics","");
    geometryType->addAttDef("Material", "Material Name", "Physics","");
    geometryType->addAttDef("Density", "Material Density", "Physics","");
    geometryType->addAttDef("State", "Material State", "Physics","");
    geometryType->addAttDef("Radlen", "Material Radiation Length", "Physics","");
    typeTree->addType(geometryType);

    eventType = heprepFactory->createHepRepType(NULL, "Event");
    eventType->addAttValue("Layer", (char*)"Event");
    typeTree->addType(eventType);

    trackType = heprepFactory->createHepRepType(eventType, "Track");
    trackType->addAttValue("Layer", (char*)"Track");

    calHitType = heprepFactory->createHepRepType(eventType, "Calorimeter Hit");
    calHitType->addAttValue("Layer", (char*)"CalHit");

    hitType = heprepFactory->createHepRepType(eventType, "Hit");
    hitType->addAttValue("Layer", (char*)"Hit");

    HepRepInstanceTree* instanceTree = heprepFactory->createHepRepInstanceTree("G4Data", "1.0", typeTree);
    heprep->addInstanceTree(instanceTree);

    parent = heprepFactory->createHepRepInstance(instanceTree, eventType);
    instanceTree->addInstance(parent);

    delete treeID;
    delete typeTree;
    delete instanceTree;
}

void G4HepRepSceneHandler::close() {
    if (writer != NULL) {
#ifdef DEBUG
        G4cout << "G4HepRepSceneHandler::close() " << G4endl;
#endif
        writer->close();
        delete writer;
        delete out;
    }
    writer = NULL;
    out = NULL;

    delete heprep;
    heprep = NULL;

    delete geometryType;
    geometryType = NULL;

    delete eventType;
    eventType = NULL;

    delete trackType;
    trackType = NULL;

    delete hitType;
    hitType = NULL;

    delete calHitType;
    calHitType = NULL;

    delete parent;
    parent = NULL;
}


void G4HepRepSceneHandler::EstablishSpecials(G4PhysicalVolumeModel& model) {
    model.DefinePointersToWorkingSpace(&currentDepth, &currentPV, &currentLV);
}

void G4HepRepSceneHandler::BeginModeling() {
    G4VSceneHandler::BeginModeling();
}

void G4HepRepSceneHandler::EndModeling() {
    G4VSceneHandler::EndModeling();
}


void G4HepRepSceneHandler::AddPrimitive (const G4Polyline& line) {

#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddPrimitive(G4Polyline&) " << line.size() << G4endl;
#endif

    HepRepFactory* factory = GetHepRepFactory();
    HepRepInstance* instance = CreateInstance(parent, trackType);

    instance->addAttValue("DrawAs", string("Line"));
    SetColour(instance, GetColour(line));
    for (size_t i=0; i < line.size(); i++) {
        G4Point3D vertex = transform * line[i];
        HepRepPoint* point = factory->createHepRepPoint(instance, vertex.x(), vertex.y(), vertex.z());
        delete point;
    }
    delete instance;
}


void G4HepRepSceneHandler::AddPrimitive (const G4Polymarker& line) {

#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddPrimitive(G4Polymarker&) " << line.size() << G4endl;
#endif
    HepRepFactory* factory = GetHepRepFactory();
    HepRepInstance* instance = CreateInstance(parent, hitType);

    instance->addAttValue("DrawAs", string("Point"));

    instance->addAttValue("MarkName", string("Square"));
    instance->addAttValue("MarkSize", 4);

    SetColour(instance, GetColour(line));
    for (size_t i=0; i < line.size(); i++) {
        G4Point3D vertex = transform * line[i];
        HepRepPoint* point = factory->createHepRepPoint(instance, vertex.x(), vertex.y(), vertex.z());
        delete point;
    }
    delete instance;
}

void G4HepRepSceneHandler::AddPrimitive (const G4Circle& circle) {
    G4Point3D center = transform * circle.GetPosition();
    G4double  radius = circle.GetWorldSize();

    HepRepFactory* factory = GetHepRepFactory();
    HepRepInstance* instance = CreateInstance(parent, hitType);

    SetColour (instance, GetColour(circle));
    instance->addAttValue("DrawAs", string("Point"));
    instance->addAttValue("MarkName", string("Dot"));
    instance->addAttValue("MarkSize", radius);

    HepRepPoint* point = factory->createHepRepPoint(instance, center.x(), center.y(), center.z());
    delete point;
    delete instance;

#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddPrimitive(G4Circle&) : center : "
        << " x : " << center.x()
        << " y : " << center.y()
        << " z : " << center.z()
        << " ,radius : " << radius << G4endl;
#endif
}



void G4HepRepSceneHandler::AddPrimitive (const G4Polyhedron& polyhedron) {

#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddPrimitive(G4Polyhedron&) " << G4endl;
#endif
    G4Normal3D surfaceNormal;
    G4Point3D vertex;

    HepRepFactory* factory = GetHepRepFactory();
    HepRepInstance* instance = CreateInstance(parent, calHitType);

    if (polyhedron.GetNoFacets()==0) return;

    G4bool notLastFace;
    do {
        HepRepInstance* face = CreateInstance(parent, calHitType);
        face->addAttValue("DrawAs", string("Polygon"));
        SetColour(face, GetColour(polyhedron));

        notLastFace = polyhedron.GetNextNormal (surfaceNormal);

        G4int edgeFlag = 1;
        G4bool notLastEdge;
        do {
            notLastEdge = polyhedron.GetNextVertex (vertex, edgeFlag);
            vertex = transform * vertex;
            HepRepPoint* point = factory->createHepRepPoint(face, vertex.x(), vertex.y(), vertex.z());
            delete point;
        } while (notLastEdge);

        delete face;

    } while (notLastFace);

    delete instance;
}


void G4HepRepSceneHandler::AddPrimitive (const G4Text& text) {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddPrimitive(G4Text&) " << G4endl;
#endif
    G4cout << "G4HepRepSceneHandler::AddPrimitive G4Text : not yet implemented. " << G4endl;
}


void G4HepRepSceneHandler::AddPrimitive (const G4Square& square) {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddPrimitive(G4Square&) " << G4endl;
#endif
    G4Point3D center = transform * square.GetPosition();
    G4double size = square.GetWorldSize();

    HepRepFactory* factory = GetHepRepFactory();
    HepRepInstance* instance = CreateInstance(parent, hitType);

    SetColour (instance, GetColour(square));
    instance->addAttValue("DrawAs", string("Point"));
    instance->addAttValue("MarkName", string("Square"));
    instance->addAttValue("MarkSize", size);

    HepRepPoint* point = factory->createHepRepPoint(instance, center.x(), center.y(), center.z());
    delete point;
    delete instance;
}


//Method for handling G4NURBS objects for drawing solids.
//Knots and Ctrl Pnts MUST be arrays of GLfloats.
void G4HepRepSceneHandler::AddPrimitive (const G4NURBS& nurb) {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddPrimitive(G4NURBS&) " << G4endl;
#endif
    G4cout << "G4HepRepSceneHandler::AddPrimitive G4NURBS : not yet implemented. " << G4endl;
}


void G4HepRepSceneHandler::AddThis (const G4Box& box) {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddThis(G4Box&) " << G4endl;
#endif
    G4VSceneHandler::AddThis (box);
}


void G4HepRepSceneHandler::AddThis (const G4Cons& cons) {
#ifdef DEBUG
  G4cout << "G4HepRepSceneHandler::AddThis(G4Cons&) " << G4endl;
#endif
  G4VSceneHandler::AddThis (cons);
}


void G4HepRepSceneHandler::AddThis (const G4Tubs& tubs) {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddThis(G4Tubs&) " << G4endl;
#endif
    G4VSceneHandler::AddThis (tubs);
}


void G4HepRepSceneHandler::AddThis (const G4Trd& trd) {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddThis(G4Trd&) " << G4endl;
#endif
    G4VSceneHandler::AddThis (trd);
}


void G4HepRepSceneHandler::AddThis (const G4Trap& trap) {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddThis(G4Trap&) " << G4endl;
#endif
    G4VSceneHandler::AddThis (trap);
}


void G4HepRepSceneHandler::AddThis (const G4Sphere& sphere) {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddThis(G4Sphere&) " << G4endl;
#endif
    G4VSceneHandler::AddThis (sphere);
}


void G4HepRepSceneHandler::AddThis (const G4Para& para) {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddThis(G4Para&) " << G4endl;
#endif
    G4VSceneHandler::AddThis (para);
}


void G4HepRepSceneHandler::AddThis (const G4Torus& torus) {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddThis(G4Torus&) " << G4endl;
#endif
    G4VSceneHandler::AddThis (torus);
}


void G4HepRepSceneHandler::AddThis (const G4Polycone& polycone) {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddThis(G4Polycone&) " << G4endl;
#endif
    G4VSceneHandler::AddThis (polycone);
}


void G4HepRepSceneHandler::AddThis (const G4Polyhedra& polyhedra) {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddThis(G4Polyhedra&) " << G4endl;
#endif
    G4VSceneHandler::AddThis (polyhedra);
}


void G4HepRepSceneHandler::AddThis (const G4VSolid& solid) {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddThis(G4VSolid&) " << G4endl;
#endif
    G4VSceneHandler::AddThis (solid);
}

void G4HepRepSceneHandler::AddThis (const G4VTrajectory& traj) {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddThis(G4VTrajectory&) " << G4endl;
#endif
    G4VSceneHandler::AddThis (traj);
}

void G4HepRepSceneHandler::AddThis (const G4VHit& hit) {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddThis(G4VHit&) " << G4endl;
#endif
    G4VSceneHandler::AddThis (hit);
}

void 
G4HepRepSceneHandler::PreAddThis 
(const G4Transform3D& objectTransformation,
			                           const G4VisAttributes& visAttribs) {

    G4VSceneHandler::PreAddThis (objectTransformation, visAttribs);

    transform = objectTransformation;
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::PreAddThis : " << objectTransformation << G4endl;
    G4cout << "G4HepRepSceneHandler::PreAddThis : " << visAttribs << G4endl;
#endif
}

void G4HepRepSceneHandler::PostAddThis () {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::PostAddThis" << G4endl;
#endif
    G4VSceneHandler::PostAddThis();
}


void G4HepRepSceneHandler::BeginPrimitives (const G4Transform3D& objectTransformation) {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::BeginPrimitives: " << G4endl;
#endif

    G4VSceneHandler::BeginPrimitives (objectTransformation);
    transform = objectTransformation;
}


void G4HepRepSceneHandler::EndPrimitives () {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::EndPrimitives" << G4endl;
#endif
    G4VSceneHandler::EndPrimitives ();
}

void G4HepRepSceneHandler::SetColour (HepRepAttribute *attribute, const G4Colour& color) {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::SetColour : red : " << color.GetRed ()   <<
                                  " green : " << color.GetGreen () <<
                                  " blue : " << color.GetBlue ()   << G4endl;
#endif
    vector<double> c;
    c.push_back(color.GetRed());
    c.push_back(color.GetGreen());
    c.push_back(color.GetBlue());
    c.push_back(color.GetAlpha());
    attribute->addAttValue("Color", c);
}

HepRepInstance* G4HepRepSceneHandler::CreateInstance(HepRepInstance* p, HepRepType* altType) {
    HepRepFactory* factory = GetHepRepFactory();
    bool event = IsEventData();
    HepRepType* type = (event) ? altType : geometryType;
    HepRepInstance* instance = factory->createHepRepInstance(p, type);
    if (!event) {
        instance->addAttValue("LVol",     currentLV->GetName());
        instance->addAttValue("Solid",    currentLV->GetSolid()->GetName());
        instance->addAttValue("EType",    currentLV->GetSolid()->GetEntityType());
        instance->addAttValue("Material", currentLV->GetMaterial()->GetName());
        instance->addAttValue("Density",  currentLV->GetMaterial()->GetDensity());
        instance->addAttValue("State",    currentLV->GetMaterial()->GetState());
        instance->addAttValue("Radlen",   currentLV->GetMaterial()->GetRadlen());
    }
    return instance;
}

bool G4HepRepSceneHandler::IsEventData () {
    return !currentPV || fReadyForTransients;
}
