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

/**
 * @author Mark Donszelmann, Joseph Perl
 */

#include <stdio.h>

#include "globals.hh"
#include <vector>
#include <strstream>
#include <iostream>
#include <iomanip>
#include <fstream>

//HepRep
#include "HEPREP/HepRep.h"

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
#include "G4VisAttributes.hh"
#include "G4VSolid.hh"
#include "G4VTrajectory.hh"
#include "G4VTrajectoryPoint.hh"
#include "G4VHit.hh"
#include "G4Scene.hh"
#include "G4Material.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"

// Streamer
#include "XMLHepRepStreamerFactory.h"

// This
#include "G4HepRep.hh"
#include "G4HepRepSceneHandler.hh"


using namespace HEPREP;

G4int G4HepRepSceneHandler::sceneCount = 0;
G4int G4HepRepSceneHandler::sceneIdCount = 0;

G4HepRepSceneHandler::G4HepRepSceneHandler (G4VGraphicsSystem& system, const G4String& name)
        : G4VSceneHandler (system, sceneIdCount++, name),
          currentDepth  (0),
          currentPV     (0),
          currentLV     (0) {

    G4cout << "G4HepRepSceneHandler::G4HepRepSceneHandler: "
           << system << " "
           << GetScene()->GetName() << G4endl;

    factory = new XMLHepRepStreamerFactory();
    writer = NULL;

    geomTypeHepRepFactory = new XMLHepRepStreamerFactory(true);
    geomInstanceHepRepFactory = new XMLHepRepStreamerFactory(true);
    geomTypeWriter = NULL;
    geomInstanceWriter = NULL;

    eventTypeHepRepFactory = new XMLHepRepStreamerFactory(true);
    eventInstanceHepRepFactory = new XMLHepRepStreamerFactory(true);
    eventTypeWriter = NULL;
    eventInstanceWriter = NULL;

    fileNo = 0;
}


G4HepRepSceneHandler::~G4HepRepSceneHandler () {
    close();

    delete factory;
    factory = NULL;

    delete geomTypeHepRepFactory;
    geomTypeHepRepFactory = NULL;
    delete geomInstanceHepRepFactory;
    geomInstanceHepRepFactory = NULL;

    delete eventTypeHepRepFactory;
    eventTypeHepRepFactory = NULL;
    delete eventInstanceHepRepFactory;
    eventInstanceHepRepFactory = NULL;
}

void G4HepRepSceneHandler::open() {
    if (writer != NULL) return;

    if (strcmp(GetScene()->GetName(), "stdout") == 0) {
        G4cout << "G4HepRepSceneHandler::open() stdout" << G4endl;
        writer = dynamic_cast<XMLHepRepStreamer*>(factory->createHepRepWriter(&G4cout));
        out  = NULL;
    } else if (strcmp(GetScene()->GetName(), "stderr") == 0) {
        G4cout << "G4HepRepSceneHandler::open() stderr" << G4endl;
        writer = dynamic_cast<XMLHepRepStreamer*>(factory->createHepRepWriter(&G4cerr));
        out = NULL;
    } else {
        G4cout << "G4HepRepSceneHandler::open() " << GetScene()->GetName() << "-" << fileNo << G4endl;
        char fileName [256];
        std::ostrstream ost(fileName, 256);
        ost << GetScene()->GetName()<< "-" << fileNo << ".heprep" << std::ends;
        out = new std::ofstream(fileName);
        writer = dynamic_cast<XMLHepRepStreamer*>(factory->createHepRepWriter(out));
    }
    fileNo++;

    // Create the HepRep that holds the Trees.
    heprep = factory->createHepRep();

    // Give this HepRep all of the layer order info for both geometry and event,
    // since these will both end up in a single HepRep.
    heprep->addLayer("Detector");
    heprep->addLayer("Event");
    heprep->addLayer("CalHit");
    heprep->addLayer("Trajectory");
    heprep->addLayer("TrajectoryPoint");
    heprep->addLayer("Hit");

    writer->write(heprep->getLayerOrder());
}


void G4HepRepSceneHandler::openGeomHepRep() {
//    G4cout << "GeomTypeWriter " << geomTypeWriter << G4endl;
//    G4cout << "GeomInstanceWriter " << geomInstanceWriter << G4endl;

    open();
    std::ostrstream geomTypeOst(geomTypeFname, 256);
    geomTypeOst << GetScene()->GetName()<< "-GeomType-" << fileNo-1 << ".tree" << std::ends;
    geomTypeOut = new std::ofstream(geomTypeFname);
    geomTypeWriter = geomTypeHepRepFactory->createHepRepWriter(geomTypeOut);

    // geometry instance writer
    std::ostrstream geomInstanceOst(geomInstanceFname, 256);
    geomInstanceOst << GetScene()->GetName() << "-GeomInstances-" << fileNo-1 << ".tree" << std::ends;
    geomInstanceOut = new std::ofstream(geomInstanceFname);
    geomInstanceWriter = geomInstanceHepRepFactory->createHepRepWriter(geomInstanceOut);

    // Create the HepRep that holds the Geometry Tree.
    geomTypeHeprep = geomTypeHepRepFactory->createHepRep();

    // Create the Geometry TypeTree.
    HepRepTreeID* geomTreeID = geomTypeHepRepFactory->createHepRepTreeID("G4GeometryTypes", "1.0");
    HepRepTypeTree* geomTypeTree = geomTypeHepRepFactory->createHepRepTypeTree(geomTreeID);
    delete geomTreeID;
    geomTypeHeprep->addTypeTree(geomTypeTree);

    // Create the top level Geometry Type.
    HepRepType* detectorType = geomTypeHepRepFactory->createHepRepType(NULL, "Detector");
//    G4cout << "Created HepRepType for Detector" << G4endl;
    geomTypeFullNameMap[G4String("Detector")] = detectorType;
    detectorType->addAttValue("Layer", G4String("Detector"));

    // Add attdefs used by all detector types.
    detectorType->addAttDef("LVol", "Logical Volume", "Physics","");
    detectorType->addAttDef("Solid", "Solid Name", "Physics","");
    detectorType->addAttDef("EType", "Entity Type", "Physics","");
    detectorType->addAttDef("Material", "Material Name", "Physics","");
    detectorType->addAttDef("Density", "Material Density", "Physics","");
    detectorType->addAttDef("State", "Material State", "Physics","");
    detectorType->addAttDef("Radlen", "Material Radiation Length", "Physics","");

    geomTypeTree->addType(detectorType);

    geomParentTypeFullNameS.push(G4String("Detector"));

    // Create the HepRep that holds the Geometry InstanceTree.
    geomInstanceHeprep = geomInstanceHepRepFactory->createHepRep();

    // Create the Goemetry InstanceTree.
    geomInstanceTree = geomInstanceHepRepFactory->
      createHepRepInstanceTree("G4GeometryData", "1.0", geomTypeTree);
    geomInstanceHeprep->addInstanceTree(geomInstanceTree);

    // Create the top level Geometry Instance.
    HepRepInstance* detectorInstance =
        geomInstanceHepRepFactory->createHepRepInstance(geomInstanceTree, detectorType);
    geomInstanceTree->addInstance(detectorInstance);
    geomParentInstanceS.push(detectorInstance);

    //delete detectorType;
    delete geomTypeTree;

    geomParentDepth = -1;
}


void G4HepRepSceneHandler::openEventHepRep() {
//    G4cout << "EventTypeWriter " << eventTypeWriter << G4endl;
//    G4cout << "EventInstanceWriter " << eventInstanceWriter << G4endl;

    open();
    if (eventTypeWriter != NULL || eventInstanceWriter != NULL) return;

    std::ostrstream eventTypeOst(eventTypeFname, 256);
    eventTypeOst << GetScene()->GetName() << "-EventTypes-" << fileNo-1 << ".tree" << std::ends;
    eventTypeOut = new std::ofstream(eventTypeFname);
    eventTypeWriter = eventTypeHepRepFactory->createHepRepWriter(eventTypeOut);

    std::ostrstream eventInstanceOst(eventInstanceFname, 256);
    eventInstanceOst << GetScene()->GetName() << "-EventInstances-" << fileNo-1 << ".tree" << std::ends;
    eventInstanceOut = new std::ofstream(eventInstanceFname);
    eventInstanceWriter = eventInstanceHepRepFactory->createHepRepWriter(eventInstanceOut);

    // Create the HepRep that holds the Event TypeTree.
    eventTypeHeprep = eventTypeHepRepFactory->createHepRep();

    // Create the Event TypeTree.
    HepRepTreeID* eventTreeID = eventTypeHepRepFactory->createHepRepTreeID("G4EventTypes", "1.0");
    HepRepTypeTree* eventTypeTree = eventTypeHepRepFactory->createHepRepTypeTree(eventTreeID);
    delete eventTreeID;
    eventTypeHeprep->addTypeTree(eventTypeTree);

    // Create the top level Event Type.
    HepRepType* eventType = eventTypeHepRepFactory->createHepRepType(NULL, "Event");
//    G4cout << "Created HepRepType for Event" << G4endl;
    eventTypeFullNameMap[G4String("Event")] = eventType;
    eventType->addAttValue("Layer", G4String("Event"));

    eventTypeTree->addType(eventType);

    eventParentTypeFullNameS.push(G4String("Event"));

    // Create the HepRep that holds the Event InstanceTree.
    eventInstanceHeprep = eventInstanceHepRepFactory->createHepRep();

    // Create the Event InstanceTree.
    eventInstanceTree = eventInstanceHepRepFactory->
        createHepRepInstanceTree("G4EventData", "1.0", eventTypeTree);
    eventInstanceHeprep->addInstanceTree(eventInstanceTree);

    // Create the top level Event Instance.
    HepRepInstance* eventInstance =
        eventInstanceHepRepFactory->createHepRepInstance(eventInstanceTree, eventType);
    eventInstanceTree->addInstance(eventInstance);
    eventParentInstanceS.push(eventInstance);

    //delete eventType;
    delete eventTypeTree;

    eventParentDepth = -1;
}


void G4HepRepSceneHandler::close() {

    if (writer == NULL) return;

    G4cout << "G4HepRepSceneHandler::close() " << G4endl;

    if (geomTypeWriter != NULL) {
        geomTypeWriter->close();
        delete geomTypeWriter;
        geomTypeWriter = NULL;

        delete geomTypeOut;
        geomTypeOut = NULL;

        geomInstanceWriter->close();
        delete geomInstanceWriter;
        geomInstanceWriter = NULL;

        delete geomInstanceOut;
        geomInstanceOut = NULL;

        delete geomTypeHeprep;
        geomTypeHeprep = NULL;

        delete geomInstanceHeprep;
        geomInstanceHeprep = NULL;

        delete geomInstanceTree;
        geomInstanceTree = NULL;

        while (!geomParentTypeFullNameS.empty())
	        geomParentTypeFullNameS.pop();

        while (!geomParentInstanceS.empty())
	        geomParentInstanceS.pop();

        geomTypeFullNameMap.clear();
    }

    if (eventTypeWriter != NULL) {
        eventTypeWriter->close();
        delete eventTypeWriter;
        eventTypeWriter = NULL;

        eventInstanceWriter->close();
        delete eventInstanceWriter;
        eventInstanceWriter = NULL;

        delete eventTypeOut;
        eventTypeOut = NULL;

        delete eventInstanceOut;
        eventInstanceOut = NULL;

        delete eventTypeHeprep;
        eventTypeHeprep = NULL;

        delete eventInstanceHeprep;
        eventInstanceHeprep = NULL;

        delete eventInstanceTree;
        eventInstanceTree = NULL;

        while (!eventParentTypeFullNameS.empty())
	        eventParentTypeFullNameS.pop();

        while (!eventParentInstanceS.empty())
	        eventParentInstanceS.pop();

        eventTypeFullNameMap.clear();
    }

    // merge files
    mergeAndDelete(geomTypeFname);
    mergeAndDelete(eventTypeFname);
    mergeAndDelete(geomInstanceFname);
    mergeAndDelete(eventInstanceFname);

    geomTypeFname[0] = 0;
    geomInstanceFname[0] = 0;
    eventTypeFname[0] = 0;
    eventInstanceFname[0] = 0;

    writer->close();
    delete writer;
    writer = NULL;

    delete out;
    out = NULL;

    delete heprep;
    heprep = NULL;
}

void G4HepRepSceneHandler::mergeAndDelete(char* fname) {
    if (fname[0] == 0) return;

    G4cout << "Merging: " << fname << G4endl;

    std::ifstream ifile( fname, std::ios::in );
    if ( ifile.fail() ) {
        G4cerr << "Failed to open: " << fname << ", will not be merged." << G4endl;
        return;
    }

    G4String line;
    while(!std::getline(ifile, line).eof()) {
        writer->printPlain(line);
        writer->printPlain("\n");
    }
    ifile.close();

    // delete file
    if (remove(fname) != 0) {
        G4cerr << "Failed to delete: " << fname << G4endl;
    }
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

    HepRepInstance* instance = CreateEventInstance(G4String("Trajectory"), 0);

    SetColour(instance, GetColour(line));

    SetLine(instance, line);

    for (size_t i=0; i < line.size(); i++) {
        G4Point3D vertex = transform * line[i];
        HepRepPoint* point =
	        eventInstanceHepRepFactory->createHepRepPoint(instance, vertex.x(), vertex.y(), vertex.z());
        delete point;
    }
    delete instance;
}


void G4HepRepSceneHandler::AddPrimitive (const G4Polymarker& line) {

#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddPrimitive(G4Polymarker&) " << line.size() << G4endl;
#endif

    HepRepInstance* instance = CreateEventInstance(G4String("Hit"), 0);

    SetColour(instance, GetColour(line));

    SetMarker(instance, line);

    // Default MarkName is set to Circle for this Type.
    switch (line.GetMarkerType()) {
        case line.dots:
            instance->addAttValue("Fill", true);
            SetColour(instance, GetColour(line), G4String("FillColor"));
            break;
        case line.circles:
            break;
        case line.squares:
            instance->addAttValue("MarkName", G4String("Box"));
            break;
        case line.line:
        default:
            instance->addAttValue("MarkName", G4String("Plus"));
            break;
    }

    for (size_t i=0; i < line.size(); i++) {
        G4Point3D vertex = transform * line[i];
        HepRepPoint* point =
	  eventInstanceHepRepFactory->createHepRepPoint(instance, vertex.x(), vertex.y(), vertex.z());
        delete point;
    }
    delete instance;
}


void G4HepRepSceneHandler::AddPrimitive (const G4Circle& circle) {

    HepRepInstance* instance = CreateEventInstance(G4String("Hit"), 0);

    G4Point3D center = transform * circle.GetPosition();

    SetColour (instance, GetColour(circle));

    SetMarker(instance, circle);

    HepRepPoint* point =
        eventInstanceHepRepFactory->createHepRepPoint(instance, center.x(), center.y(), center.z());
    delete point;
    delete instance;

#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddPrimitive(G4Circle&) : center : "
        << " x : " << center.x()
        << " y : " << center.y()
        << " z : " << center.z()
        << " ,radius : " << size << G4endl;
#endif
}


void G4HepRepSceneHandler::AddPrimitive (const G4Polyhedron& polyhedron) {

#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddPrimitive(G4Polyhedron&) " << G4endl;
#endif
    G4Normal3D surfaceNormal;
    G4Point3D vertex;

    if (polyhedron.GetNoFacets()==0) return;

//    G4cout << "Will use ParentDepth of " << geomParentDepth << G4endl;

    HepRepInstance* instance;
    if (IsEventData())
        instance = CreateEventInstance(G4String("CalHit"), 0);
    else
        instance = CreateGeomInstance(currentLV->GetName(), currentDepth);

    G4bool notLastFace;
    do {
        HepRepInstance* face;
        if (IsEventData())
	        face = CreateEventInstance(G4String("Face"), 1);
        else
	        face = CreateGeomInstance(G4String("Face"), currentDepth+1);

        SetLine(face, polyhedron);
        SetColour(face, GetColour(polyhedron));
        if (IsEventData()) SetColour(face, GetColour(polyhedron), G4String("FillColor"));

        notLastFace = polyhedron.GetNextNormal (surfaceNormal);

        G4int edgeFlag = 1;
        G4bool notLastEdge;
        do {
	        notLastEdge = polyhedron.GetNextVertex (vertex, edgeFlag);
            vertex = transform * vertex;
	        HepRepPoint* point;
	        if (IsEventData())
	            point = eventInstanceHepRepFactory->createHepRepPoint(face, vertex.x(), vertex.y(), vertex.z());
	        else
	            point = geomInstanceHepRepFactory->createHepRepPoint(face, vertex.x(), vertex.y(), vertex.z());
	        delete point;
        } while (notLastEdge);

        delete face;

    } while (notLastFace);
    delete instance;
}


void G4HepRepSceneHandler::AddPrimitive (const G4Text&) {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddPrimitive(G4Text&) " << G4endl;
#endif
    G4cout << "G4HepRepSceneHandler::AddPrimitive G4Text : not yet implemented. " << G4endl;
}


void G4HepRepSceneHandler::AddPrimitive (const G4Square& square) {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddPrimitive(G4Square&) " << G4endl;
#endif

    HepRepInstance* instance = CreateEventInstance(G4String("Hit"), 0);

    G4Point3D center = transform * square.GetPosition();

    SetColour (instance, GetColour(square));

    SetMarker(instance, square);

    instance->addAttValue("MarkName", G4String("Box"));

    HepRepPoint* point =
        eventInstanceHepRepFactory->createHepRepPoint(instance, center.x(), center.y(), center.z());
    delete point;
    delete instance;
}


//Method for handling G4NURBS objects for drawing solids.
//Knots and Ctrl Pnts MUST be arrays of GLfloats.
void G4HepRepSceneHandler::AddPrimitive (const G4NURBS&) {
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

    std::vector<G4AttValue>* trajAttValues = traj.CreateAttValues();
    const std::map<G4String,G4AttDef>* trajAttDefs = traj.GetAttDefs();

    HepRepInstance* trajInstance = CreateEventInstance(G4String("Trajectory"), 0, trajAttDefs, trajAttValues);

    // Set the LineColor attribute according to the particle charge.
    float redness = 0.;
    float greenness = 0.;
    float blueness = 0.;
    const G4double charge = traj.GetCharge();
    if(charge>0.)      blueness =  1.; // Blue = positive.
    else if(charge<0.) redness  =  1.; // Red = negative.
    else               greenness = 1.; // Green = neutral.
    // Not sure how to set color in this system
    //trajInstance->addAttValue(new G4String("Color"), redness, greenness, blueness);

    // Specify the polyline by using the trajectory points.
    G4int i;
    for (i = 0; i < traj.GetPointEntries(); i++) {
        G4VTrajectoryPoint* aTrajectoryPoint = traj.GetPoint(i);
        G4Point3D vertex = aTrajectoryPoint->GetPosition();
        HepRepPoint* point =
            eventInstanceHepRepFactory->createHepRepPoint(trajInstance, vertex.x(), vertex.y(), vertex.z());
        delete point;
    }

    for (i = 0; i < traj.GetPointEntries(); i++) {
        G4VTrajectoryPoint* aTrajectoryPoint = traj.GetPoint(i);

        std::vector<G4AttValue>* pointAttValues = aTrajectoryPoint->CreateAttValues();
        const std::map<G4String,G4AttDef>* pointAttDefs = aTrajectoryPoint->GetAttDefs();

        HepRepInstance* pointInstance = CreateEventInstance(G4String("Point"), 1, pointAttDefs, pointAttValues);

        // Specify the position of the trajectory point.
        G4Point3D vertex = aTrajectoryPoint->GetPosition();
        HepRepPoint* point =
            eventInstanceHepRepFactory->createHepRepPoint(pointInstance, vertex.x(), vertex.y(), vertex.z());
        delete point;
        delete pointInstance;
    }
    delete trajInstance;
}


void G4HepRepSceneHandler::AddThis (const G4VHit& hit) {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::AddThis(G4VHit&) " << G4endl;
#endif
    G4VSceneHandler::AddThis (hit);
}


void
G4HepRepSceneHandler::PreAddThis (const G4Transform3D& objectTransformation,
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


void G4HepRepSceneHandler::SetColour (HepRepAttribute *attribute,
				      const G4Colour& color,
				      const G4String& key) {
#ifdef DEBUG
    G4cout << "G4HepRepSceneHandler::SetColour : red : " << color.GetRed ()   <<
                                  " green : " << color.GetGreen () <<
                                  " blue : " << color.GetBlue ()   << G4endl;
#endif

    attribute->addAttValue(key, color.GetRed(), color.GetGreen(), color.GetBlue(),color.GetAlpha());
}


void G4HepRepSceneHandler::SetLine (HepRepInstance *instance, const G4Visible& visible) {
    G4VisAttributes atts = visible.GetVisAttributes();
    instance->addAttValue("LineWidth", atts.GetLineWidth());
    switch (atts.GetLineStyle()) {
        case G4VisAttributes::dotted:
            instance->addAttValue("LineStyle", G4String("Dotted"));
            break;
        case G4VisAttributes::dashed:
            instance->addAttValue("LineStyle", G4String("Dashed"));
            break;
        case G4VisAttributes::unbroken:
        default:
            break;
    }
}


void G4HepRepSceneHandler::SetMarker (HepRepInstance *instance, const G4VMarker& marker) {
    MarkerSizeType markerType;
    G4double size = GetMarkerRadius( marker , markerType );
    instance->addAttValue("MarkSize", size);

    if (markerType == screen)
	    instance->addAttValue("MarkType", G4String("Symbol"));
    if (marker.GetFillStyle() == G4VMarker::noFill) {
        instance->addAttValue("Fill", false);
    } else {
        SetColour(instance, GetColour(marker), G4String("FillColor"));
    }
}


HepRepInstance* G4HepRepSceneHandler::CreateGeomInstance(G4String typeName, G4int depth) {

    if (geomInstanceWriter == NULL)
        openGeomHepRep();

//    G4cout << "typeName " << typeName << G4endl;
//    G4cout << "depth " << depth << G4endl;
//    G4cout << "parent Depth " << geomParentDepth << G4endl;

    // Pop stacks to get current parent onto the top
    while (depth <= geomParentDepth-- ) {
//        G4cout << "Popping old parent of fullTypeName "
//               << geomParentTypeFullNameS.top() << G4endl;
        geomParentTypeFullNameS.pop();
        geomParentInstanceS.pop();
    }

    // Create the full type name.
    G4String typeFullName = geomParentTypeFullNameS.top() + "/" + typeName;
//    G4cout << "Constructed typeFullName as " << typeFullName << G4endl;

    // See if typeFullName is already in map
    HepRepType* type = geomTypeFullNameMap[typeFullName];

    // Add a type if this one is not yet in the map
    if (type==NULL) {
//        G4cout << "Type needs to be added to the map" << G4endl;
//        G4cout << "parentTypeFullNameS.size() = " << geomParentTypeFullNameS.size() << G4endl;
//        G4cout << "parentTypeFullNameS.top() = " << geomParentTypeFullNameS.top() << G4endl;
        HepRepType* geomParentType = geomTypeFullNameMap[geomParentTypeFullNameS.top()];
//        G4cout << "parentType was found to be " << geomParentType << G4endl;
//        G4cout << "parentType name = " << geomParentType->getName() << G4endl;
//        G4cout << "typeFactory = " << geomTypeHepRepFactory << G4endl;
//        G4cout << "Creating HepRepType for " << typeName << G4endl;
        type = geomTypeHepRepFactory->createHepRepType(geomParentType, typeFullName);
//        G4cout << "Created new type " << type << G4endl;
        geomTypeFullNameMap[typeFullName] = type;
//        G4cout << "Entered new type into map " << G4endl;

        // Don't need to specify any attribute definitions, since all Detector types inherit
        // the same set of definitions from the top level detector type.

        // Add type-specific attributes
        if (typeName == G4String("Face"))
        	type->addAttValue("DrawAs", G4String("Polygon"));
    }

    // add the instance
//    G4cout << "Adding instance of type " << type->getName() <<
//              " and geomParentInstance "  << geomParentInstanceS.top() << G4endl;
    HepRepInstance* instance =
      geomInstanceHepRepFactory->createHepRepInstance(geomParentInstanceS.top(), type);
//    G4cout << "Added instance " << instance << G4endl;

    // add the latest info to the stacks
    geomParentDepth = depth;
    geomParentTypeFullNameS.push(typeFullName);
    geomParentInstanceS.push(instance);

    // add geometry attributes
    instance->addAttValue("LVol",     currentLV->GetName());
    instance->addAttValue("Solid",    currentLV->GetSolid()->GetName());
    instance->addAttValue("EType",    currentLV->GetSolid()->GetEntityType());
    instance->addAttValue("Material", currentLV->GetMaterial()->GetName());
    instance->addAttValue("Density",  currentLV->GetMaterial()->GetDensity());
    instance->addAttValue("Radlen",   currentLV->GetMaterial()->GetRadlen());

    switch (currentLV->GetMaterial()->GetState()) {
    case kStateSolid:
        instance->addAttValue("State", G4String("Solid"));
        break;
    case kStateLiquid:
        instance->addAttValue("State", G4String("Liquid"));
        break;
    case kStateGas:
        instance->addAttValue("State", G4String("Gas"));
        break;
    case kStateUndefined:
    default:
        instance->addAttValue("State", G4String("Undefined"));
        break;
    }
    return instance;
}


HepRepInstance* G4HepRepSceneHandler::CreateEventInstance(G4String typeName, G4int depth,
                                                          const std::map<G4String,G4AttDef>* attDefs,
                                                          std::vector<G4AttValue>* attValues) {

    if (eventInstanceWriter == NULL)
        openEventHepRep();

//    G4cout << "typeName " << typeName << G4endl;
//    G4cout << "depth " << depth << G4endl;
//    G4cout << "parent Depth " << eventParentDepth << G4endl;

    // Pop stacks to get current parent onto the top
    while (depth <= eventParentDepth-- ) {
//        G4cout << "Popping old parent of fullTypeName "
//               << eventParentTypeFullNameS.top() << G4endl;
        eventParentTypeFullNameS.pop();
        eventParentInstanceS.pop();
    }

    // Create the full type name.
    G4String typeFullName = eventParentTypeFullNameS.top() + "/" + typeName;
//    G4cout << "Constructed typeFullName as " << typeFullName << G4endl;

    // See if typeFullName is already in map
    HepRepType* type = eventTypeFullNameMap[typeFullName];

    // Add a type if this one is not yet in the map
    if (type==NULL) {
//        G4cout << "Type needs to be added to the map" << G4endl;
//        G4cout << "parentTypeFullNameS.size() = " << eventParentTypeFullNameS.size() << G4endl;
//        G4cout << "parentTypeFullNameS.top() = " << eventParentTypeFullNameS.top() << G4endl;
        HepRepType* eventParentType = eventTypeFullNameMap[eventParentTypeFullNameS.top()];
//        G4cout << "parentType was found to be " << eventParentType << G4endl;
//        G4cout << "parentType name = " << eventParentType->getName() << G4endl;
//        G4cout << "typeFactory = " << eventTypeHepRepFactory << G4endl;
//        G4cout << "Creating HepRepType for " << typeName << G4endl;
        type = eventTypeHepRepFactory->createHepRepType(eventParentType, typeFullName);
//        G4cout << "Created new type " << type << G4endl;
        eventTypeFullNameMap[typeFullName] = type;
//        G4cout << "Entered new type into map " << G4endl;

        // Specify additional attribute definitions.  Assumes that all attributes for
        //  the given type can be found from the first instance of this type.
        std::vector<G4AttValue>::iterator iAttVal;
        std::map<G4String,G4AttDef>::const_iterator iAttDef;
        if (attValues && attDefs) {
	        for (iAttVal = attValues->begin(); iAttVal != attValues->end(); ++iAttVal) {
	            iAttDef = attDefs->find(iAttVal->GetName());
	            if (iAttDef != attDefs->end()) {
            	    // Protect against incorrect use of Category.  Anything value other than the
            	    // standard ones will be considered to be in the physics category.
            	    G4String category = iAttDef->second.GetCategory();
            	    if (strcmp(category,"Draw")!=0 &&
                		strcmp(category,"Physics")!=0 &&
                		strcmp(category,"Association")!=0 &&
                		strcmp(category,"PickAction")!=0) {

	                    category = "Physics";
	                }
	                type->addAttDef(iAttVal->GetName(), iAttDef->second.GetDesc(),
			        category, iAttDef->second.GetExtra());
	            }
	        }
        }

        // Add type-specific attributes
        if (typeName == G4String("CalHit")) {
        	type->addAttValue("Layer", G4String("CalHit"));
        	type->addAttValue("Fill", true);
        } else if (typeName == G4String("Trajectory")) {
        	type->addAttValue("Layer", G4String("Trajectory"));
        	type->addAttValue("DrawAs", G4String("Line"));
        } else if (typeName == G4String("Point")) {
        	type->addAttValue("Layer", G4String("TrajectoryPoint"));
        	type->addAttValue("MarkName", G4String("Circle"));
        	type->addAttValue("MarkType", G4String("Real"));
        	type->addAttValue("Fill", true);
        	type->addAttValue("Visibility",true);
        } else if (typeName == G4String("Hit")) {
        	type->addAttValue("Layer", G4String("Hit"));
        	type->addAttValue("DrawAs", G4String("Point"));
        	type->addAttValue("MarkName", G4String("Circle"));
        	type->addAttValue("MarkType", G4String("Real"));
        	type->addAttValue("Fill", true);
        } else if (typeName == G4String("Face")) {
        	type->addAttValue("Layer", G4String("CalHit"));
        	type->addAttValue("DrawAs", G4String("Line"));
        }
    }

    // add the instance
//    G4cout << "Adding instance of type " << type->getName() <<
//              " and eventParentInstance "  << eventParentInstanceS.top() << G4endl;
    HepRepInstance* instance =
        eventInstanceHepRepFactory->createHepRepInstance(eventParentInstanceS.top(), type);
//    G4cout << "Added instance " << instance << G4endl;

    // add the latest info to the stacks
    eventParentDepth = depth;
    eventParentTypeFullNameS.push(typeFullName);
    eventParentInstanceS.push(instance);

    // Copy the instance's G4AttValues to HepRepAttValues.
    std::vector<G4AttValue>::iterator iAttVal;
    if (attValues && attDefs) {
        for (iAttVal = attValues->begin(); iAttVal != attValues->end(); ++iAttVal) {
	        std::map<G4String,G4AttDef>::const_iterator iAttDef =
	        attDefs->find(iAttVal->GetName());
	        if (iAttDef == attDefs->end()) {
//	            G4cout << "G4HepRepFileSceneHandler::AddThis(traj):"
//                          "\n  WARNING: no matching definition for attribute \""
//                       << iAttVal->GetName() << "\", value: "
//                       << iAttVal->GetValue();
	        } else {
                // Use GetDesc rather than GetName once WIRED can handle names with spaces in them.
                //instance->addAttValue(iAttDef->second.GetDesc(), iAttVal->GetValue());
                instance->addAttValue(iAttVal->GetName(), iAttVal->GetValue());
        	}
        }
        delete attValues;  // AttValues must be deleted after use.
    }

    return instance;
}

bool G4HepRepSceneHandler::IsEventData () {
    return !currentPV || fReadyForTransients;
}
