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
#include <iostream>
#include <sstream>
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

// JHepRep
#include "XMLHepRepFactory.h"

// This
#include "G4HepRep.hh"
#include "G4HepRepSceneHandler.hh"
#include "G4HepRepViewer.hh"


using namespace HEPREP;
using namespace std;

G4int G4HepRepSceneHandler::sceneCount = 0;

G4HepRepSceneHandler::G4HepRepSceneHandler (G4VGraphicsSystem& system, const G4String& name)
        : G4VSceneHandler (system, sceneCount++, name),
          eventNumber   (1),
          currentDepth  (0),
          currentPV     (0),
          currentLV     (0),
          heprep	(NULL) {

#ifdef DEBUG
    cout << "G4HepRepSceneHandler::G4HepRepSceneHandler: " << system << endl;
//           << GetScene()->GetName() << endl;
#endif

    factory = new XMLHepRepFactory();
    writer = NULL;

    Open(GetScene() == NULL ? G4String("G4HepRepOutput.zip") : GetScene()->GetName());
    OpenHepRep();
}


G4HepRepSceneHandler::~G4HepRepSceneHandler () {
#ifdef DEBUG
    cout << "G4HepRepSceneHandler::~G4HepRepSceneHandler() " << endl;
#endif
    Close();

    delete factory;
    factory = NULL;

    dynamic_cast<G4HepRep*>(GetGraphicsSystem())->RemoveSceneHandler();
}

void G4HepRepSceneHandler::Open(G4String name) {
    if (writer != NULL) return;

    if (strcmp(name, "stdout") == 0) {
#ifdef DEBUG
        cout << "G4HepRepSceneHandler::Open() stdout" << endl;
#endif
        writer = factory->createHepRepWriter(&cout, false, false);
        out  = NULL;
    } else if (strcmp(name, "stderr") == 0) {
#ifdef DEBUG
        cout << "G4HepRepSceneHandler::Open() stderr" << endl;
#endif
        writer = factory->createHepRepWriter(&cerr, false, false);
        out = NULL;
    } else {
#ifdef DEBUG
        cout << "G4HepRepSceneHandler::Open() " << name << endl;
#endif
        out = new ofstream(name.c_str(), std::ios::out | std::ios::binary );
        bool zip = name.substr(name.size()-4, 4) == G4String(".zip");
        bool gz = name.substr(name.size()-3, 3) == G4String(".gz");
        writer = factory->createHepRepWriter(out, zip, zip || gz);
    }
}


void G4HepRepSceneHandler::OpenHepRep() {
#ifdef DEBUG
    cout << "G4HepRepSceneHandler::OpenHepRep() " << endl;
#endif

    if (heprep != NULL) return;
//    cout << "GeomTypeWriter " << geomTypeWriter << endl;
//    cout << "GeomInstanceWriter " << geomInstanceWriter << endl;

    heprepEmpty = true;

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

    // Create the Geometry TypeTree.
    HepRepTreeID* geometryTreeID = factory->createHepRepTreeID("G4GeometryTypes", "1.0");
    HepRepTypeTree* geometryTypeTree = factory->createHepRepTypeTree(geometryTreeID);
    heprep->addTypeTree(geometryTypeTree);

    // Create the top level Geometry Type.
    detectorType = factory->createHepRepType(NULL, "Detector");
//    cout << "Created HepRepType for Detector" << endl;
//    geomTypeFullNameMap[G4String("Detector")] = detectorType;
    detectorType->addAttValue("Layer", G4String("Detector"));

    // Add attdefs used by all detector types.
    detectorType->addAttDef("LVol", "Logical Volume", "Physics","");
    detectorType->addAttDef("Solid", "Solid Name", "Physics","");
    detectorType->addAttDef("EType", "Entity Type", "Physics","");
    detectorType->addAttDef("Material", "Material Name", "Physics","");
    detectorType->addAttDef("Density", "Material Density", "Physics","");
    detectorType->addAttDef("State", "Material State", "Physics","");
    detectorType->addAttDef("Radlen", "Material Radiation Length", "Physics","");

    geometryTypeTree->addType(detectorType);

//    geomParentTypeFullNameS.push(G4String("Detector"));

    // Create the Goemetry InstanceTree.
    geometryInstanceTree = factory->createHepRepInstanceTree("G4GeometryData", "1.0", geometryTypeTree);
    heprep->addInstanceTree(geometryInstanceTree);

    // Create the top level Geometry Instance.
    HepRepInstance* detectorInstance = factory->createHepRepInstance(geometryInstanceTree, detectorType);
//    geomParentInstanceS.push(detectorInstance);

//    geomParentDepth = -1;

    // Create the Event TypeTree.
    HepRepTreeID* eventTreeID = factory->createHepRepTreeID("G4EventTypes", "1.0");
    HepRepTypeTree* eventTypeTree = factory->createHepRepTypeTree(eventTreeID);
    heprep->addTypeTree(eventTypeTree);

    // Create the top level Event Type.
    eventType = factory->createHepRepType(NULL, "Event");
//    cout << "Created HepRepType for Event" << endl;
//    eventTypeFullNameMap[G4String("Event")] = eventType;
    eventType->addAttValue("Layer", G4String("Event"));

    eventTypeTree->addType(eventType);

//    eventParentTypeFullNameS.push(G4String("Event"));

    // Create the Event InstanceTree.
    eventInstanceTree = factory->createHepRepInstanceTree("G4EventData", "1.0", eventTypeTree);
    heprep->addInstanceTree(eventInstanceTree);

    // Create the top level Event Instance.
    HepRepInstance* eventInstance = factory->createHepRepInstance(eventInstanceTree, eventType);
//    eventParentInstanceS.push(eventInstance);

//    eventParentDepth = -1;

}

/**
 * Returns true if the HepRep was (already) closed, false if the HepRep is still open
 */
bool G4HepRepSceneHandler::CloseHepRep() {

#ifdef DEBUG
    cout << "G4HepRepSceneHandler::CloseHepRep() " << endl;
#endif

    if (heprep == NULL) return true;

    if (!heprepEmpty) {
        // add geometry to the heprep
        G4HepRepViewer* viewer = dynamic_cast<G4HepRepViewer*>(GetCurrentViewer());
        viewer->ProcessScene();

        // write out the heprep
        char eventName[255];
        sprintf(eventName, "event%010d.heprep", eventNumber);
        writer->write(heprep, eventName);
        eventNumber++;
    }

    delete heprep;
    heprep = NULL;
    heprepEmpty = true;

    return true;
}

void G4HepRepSceneHandler::Close() {

#ifdef DEBUG
    cout << "G4HepRepSceneHandler::Close() " << endl;
#endif

    if (writer == NULL) return;

    CloseHepRep();

    writer->close();
    delete writer;
    writer = NULL;

    delete out;
    out = NULL;
}

void G4HepRepSceneHandler::EstablishSpecials(G4PhysicalVolumeModel& model) {
#ifdef DEBUG
    cout << "G4HepRepSceneHandler::EstablishSpecials(G4PhysicalVolumeModel&) " << endl;
#endif
    model.DefinePointersToWorkingSpace(&currentDepth, &currentPV, &currentLV);
}


void G4HepRepSceneHandler::BeginModeling() {
#ifdef DEBUG
    cout << "G4HepRepSceneHandler::BeginModeling() " << endl;
#endif
    G4VSceneHandler::BeginModeling();
}


void G4HepRepSceneHandler::EndModeling() {
#ifdef DEBUG
    cout << "G4HepRepSceneHandler::EndModeling() " << endl;
#endif
    G4VSceneHandler::EndModeling();
}


void G4HepRepSceneHandler::AddPrimitive (const G4Polyline& line) {

#ifdef PDEBUG
    cout << "G4HepRepSceneHandler::AddPrimitive(G4Polyline&) " << line.size() << endl;
#endif

    HepRepInstance* instance = CreateEventInstance(G4String("Trajectory"), 0);

    SetColour(instance, GetColour(line));

    SetLine(instance, line);

    for (size_t i=0; i < line.size(); i++) {
        G4Point3D vertex = transform * line[i];
        factory->createHepRepPoint(instance, vertex.x(), vertex.y(), vertex.z());
    }
}


void G4HepRepSceneHandler::AddPrimitive (const G4Polymarker& line) {

#ifdef PDEBUG
    cout << "G4HepRepSceneHandler::AddPrimitive(G4Polymarker&) " << line.size() << endl;
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
        factory->createHepRepPoint(instance, vertex.x(), vertex.y(), vertex.z());
    }
}


void G4HepRepSceneHandler::AddPrimitive (const G4Circle& circle) {
#ifdef PDEBUG
    cout << "G4HepRepSceneHandler::AddPrimitive(G4Circle&) " << endl;
#endif

    HepRepInstance* instance = CreateEventInstance(G4String("Hit"), 0);

    G4Point3D center = transform * circle.GetPosition();

    SetColour (instance, GetColour(circle));

    SetMarker(instance, circle);

    factory->createHepRepPoint(instance, center.x(), center.y(), center.z());
}


void G4HepRepSceneHandler::AddPrimitive (const G4Polyhedron& polyhedron) {

#ifdef PDEBUG
    cout << "G4HepRepSceneHandler::AddPrimitive(G4Polyhedron&) " << endl;
#endif
    G4Normal3D surfaceNormal;
    G4Point3D vertex;

    if (polyhedron.GetNoFacets()==0) return;

//    cout << "Will use ParentDepth of " << geomParentDepth << endl;

    HepRepInstance* instance;
    if (IsEventData()) {
        instance = CreateEventInstance(G4String("CalHit"), 0);
    } else {
        instance = CreateGeometryInstance(currentLV->GetName(), currentDepth);
    }

    G4bool notLastFace;
    do {
        HepRepInstance* face;
        if (IsEventData()) {
	        face = CreateEventInstance(G4String("Face"), 1);
        } else {
	        face = CreateGeometryInstance(G4String("Face"), currentDepth+1);
        }

        SetLine(face, polyhedron);
        SetColour(face, GetColour(polyhedron));
        if (IsEventData()) SetColour(face, GetColour(polyhedron), G4String("FillColor"));

        notLastFace = polyhedron.GetNextNormal (surfaceNormal);

        G4int edgeFlag = 1;
        G4bool notLastEdge;
        do {
	        notLastEdge = polyhedron.GetNextVertex (vertex, edgeFlag);
            vertex = transform * vertex;
	        if (IsEventData()) {
	            factory->createHepRepPoint(face, vertex.x(), vertex.y(), vertex.z());
	        } else {
	            factory->createHepRepPoint(face, vertex.x(), vertex.y(), vertex.z());
	        }
        } while (notLastEdge);
    } while (notLastFace);
}


void G4HepRepSceneHandler::AddPrimitive (const G4Text&) {
#ifdef PDEBUG
    cout << "G4HepRepSceneHandler::AddPrimitive(G4Text&) " << endl;
#endif
    cout << "G4HepRepSceneHandler::AddPrimitive G4Text : not yet implemented. " << endl;
}


void G4HepRepSceneHandler::AddPrimitive (const G4Square& square) {
#ifdef PDEBUG
    cout << "G4HepRepSceneHandler::AddPrimitive(G4Square&) " << endl;
#endif

    HepRepInstance* instance = CreateEventInstance(G4String("Hit"), 0);

    G4Point3D center = transform * square.GetPosition();

    SetColour (instance, GetColour(square));

    SetMarker(instance, square);

    instance->addAttValue("MarkName", G4String("Box"));

    factory->createHepRepPoint(instance, center.x(), center.y(), center.z());
}


//Method for handling G4NURBS objects for drawing solids.
//Knots and Ctrl Pnts MUST be arrays of GLfloats.
void G4HepRepSceneHandler::AddPrimitive (const G4NURBS&) {
#ifdef PDEBUG
    cout << "G4HepRepSceneHandler::AddPrimitive(G4NURBS&) " << endl;
#endif
    cout << "G4HepRepSceneHandler::AddPrimitive G4NURBS : not yet implemented. " << endl;
}


void G4HepRepSceneHandler::AddThis (const G4VTrajectory& traj) {
#ifdef PDEBUG
    cout << "G4HepRepSceneHandler::AddThis(G4VTrajectory&) " << endl;
#endif

    vector<G4AttValue>* trajAttValues = traj.CreateAttValues();
    const map<G4String,G4AttDef>* trajAttDefs = traj.GetAttDefs();

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
        factory->createHepRepPoint(trajInstance, vertex.x(), vertex.y(), vertex.z());
    }

    for (i = 0; i < traj.GetPointEntries(); i++) {
        G4VTrajectoryPoint* aTrajectoryPoint = traj.GetPoint(i);

        vector<G4AttValue>* pointAttValues = aTrajectoryPoint->CreateAttValues();
        const map<G4String,G4AttDef>* pointAttDefs = aTrajectoryPoint->GetAttDefs();

        HepRepInstance* pointInstance = CreateEventInstance(G4String("Point"), 1, pointAttDefs, pointAttValues);

        // created
        delete pointAttValues;

        // Specify the position of the trajectory point.
        G4Point3D vertex = aTrajectoryPoint->GetPosition();
        factory->createHepRepPoint(pointInstance, vertex.x(), vertex.y(), vertex.z());
    }
}



void
G4HepRepSceneHandler::PreAddThis (const G4Transform3D& objectTransformation,
				  const G4VisAttributes& visAttribs) {

    G4VSceneHandler::PreAddThis (objectTransformation, visAttribs);

    transform = objectTransformation;
#ifdef DEBUG
    cout << "G4HepRepSceneHandler::PreAddThis(G4Transform3D&, G4VisAttributes&)" << endl;
#endif
}


void G4HepRepSceneHandler::PostAddThis () {
#ifdef DEBUG
    cout << "G4HepRepSceneHandler::PostAddThis()" << endl;
#endif
    G4VSceneHandler::PostAddThis();
}


void G4HepRepSceneHandler::BeginPrimitives (const G4Transform3D& objectTransformation) {
#ifdef DEBUG
    cout << "G4HepRepSceneHandler::BeginPrimitives(G4Transform3D&)" << endl;
#endif

    G4VSceneHandler::BeginPrimitives (objectTransformation);
    transform = objectTransformation;
    heprepEmpty = false;
}


void G4HepRepSceneHandler::EndPrimitives () {
#ifdef DEBUG
    cout << "G4HepRepSceneHandler::EndPrimitives" << endl;
#endif
    G4VSceneHandler::EndPrimitives ();
}


void G4HepRepSceneHandler::SetColour (HepRepAttribute *attribute,
				      const G4Colour& color,
				      const G4String& key) {
#ifdef VDEBUG
    cout << "G4HepRepSceneHandler::SetColour : red : " << color.GetRed ()   <<
                                  " green : " << color.GetGreen () <<
                                  " blue : " << color.GetBlue ()   << endl;
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

    if (markerType == screen) instance->addAttValue("MarkType", G4String("Symbol"));
    if (marker.GetFillStyle() == G4VMarker::noFill) {
        instance->addAttValue("Fill", false);
    } else {
        SetColour(instance, GetColour(marker), G4String("FillColor"));
    }
}


HepRepInstance* G4HepRepSceneHandler::CreateGeometryInstance(G4String typeName, G4int depth) {

    cout << "typeName " << typeName << endl;
    cout << "depth " << depth << endl;

    return factory->createHepRepInstance(geometryInstanceTree, detectorType);
/*

    // Pop stacks to get current parent onto the top
//    while (depth <= geomParentDepth-- ) {
//        cout << "Popping old parent of fullTypeName "
//               << geomParentTypeFullNameS.top() << endl;
//        geomParentTypeFullNameS.pop();
//        geomParentInstanceS.pop();
//    }

    // Create the full type name.
    G4String typePath = geomParentTypeFullNameS.top() + "/" + typeName;
//    cout << "Constructed typeFullName as " << typeFullName << endl;

    // See if typeFullName is already in map
    HepRepType* type = geometryTypeByPath[typePath];

    // Add a type if this one is not yet in the map
    if (type==NULL) {
//        cout << "Type needs to be added to the map" << endl;
//        cout << "parentTypeFullNameS.size() = " << geomParentTypeFullNameS.size() << endl;
//        cout << "parentTypeFullNameS.top() = " << geomParentTypeFullNameS.top() << endl;
        HepRepType* geometryParentType = geometryTypeByDepth[depth];
//        cout << "parentType was found to be " << geomParentType << endl;
//        cout << "parentType name = " << geomParentType->getName() << endl;
//        cout << "typeFactory = " << geomTypeHepRepFactory << endl;
//        cout << "Creating HepRepType for " << typeName << endl;
        type = factory->createHepRepType(geometryParentType, typePath);
//        cout << "Created new type " << type << endl;
        geometryTypeByPath[typePath] = type;
//        cout << "Entered new type into map " << endl;

        // Don't need to specify any attribute definitions, since all Detector types inherit
        // the same set of definitions from the top level detector type.

        // Add type-specific attributes
        if (typeName == G4String("Face")) type->addAttValue("DrawAs", G4String("Polygon"));
    }

    // add the instance
//    cout << "Adding instance of type " << type->getName() <<
//              " and geomParentInstance "  << geomParentInstanceS.top() << endl;
    HepRepInstance* instance = factory->createHepRepInstance(geomParentInstanceS.top(), type);
//    cout << "Added instance " << instance << endl;

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
*/
}


HepRepInstance* G4HepRepSceneHandler::CreateEventInstance(G4String typeName, G4int depth,
                                                          const map<G4String,G4AttDef>* /* attDefs */,
                                                          vector<G4AttValue>* /* attValues */) {

    cout << "typeName " << typeName << endl;
    cout << "depth " << depth << endl;

    return factory->createHepRepInstance(eventInstanceTree, eventType);
/*
    // Pop stacks to get current parent onto the top
    while (depth <= eventParentDepth-- ) {
//        cout << "Popping old parent of fullTypeName "
//               << eventParentTypeFullNameS.top() << endl;
        eventParentTypeFullNameS.pop();
        eventParentInstanceS.pop();
    }

    // Create the full type name.
    G4String typeFullName = eventParentTypeFullNameS.top() + "/" + typeName;
//    cout << "Constructed typeFullName as " << typeFullName << endl;

    // See if typeFullName is already in map
    HepRepType* type = eventTypeFullNameMap[typeFullName];

    // Add a type if this one is not yet in the map
    if (type==NULL) {
//        cout << "Type needs to be added to the map" << endl;
//        cout << "parentTypeFullNameS.size() = " << eventParentTypeFullNameS.size() << endl;
//        cout << "parentTypeFullNameS.top() = " << eventParentTypeFullNameS.top() << endl;
        HepRepType* eventParentType = eventTypeFullNameMap[eventParentTypeFullNameS.top()];
//        cout << "parentType was found to be " << eventParentType << endl;
//        cout << "parentType name = " << eventParentType->getName() << endl;
//        cout << "typeFactory = " << eventTypeHepRepFactory << endl;
//        cout << "Creating HepRepType for " << typeName << endl;
        type = factory->createHepRepType(eventParentType, typeFullName);
//        cout << "Created new type " << type << endl;
        eventTypeFullNameMap[typeFullName] = type;
//        cout << "Entered new type into map " << endl;

        // Specify additional attribute definitions.  Assumes that all attributes for
        //  the given type can be found from the first instance of this type.
        vector<G4AttValue>::iterator iAttVal;
        map<G4String,G4AttDef>::const_iterator iAttDef;
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
//    cout << "Adding instance of type " << type->getName() <<
//              " and eventParentInstance "  << eventParentInstanceS.top() << endl;
    HepRepInstance* instance = factory->createHepRepInstance(eventParentInstanceS.top(), type);
//    cout << "Added instance " << instance << endl;

    // add the latest info to the stacks
    eventParentDepth = depth;
    eventParentTypeFullNameS.push(typeFullName);
    eventParentInstanceS.push(instance);

    // Copy the instance's G4AttValues to HepRepAttValues.
    vector<G4AttValue>::iterator iAttVal;
    if (attValues && attDefs) {
        for (iAttVal = attValues->begin(); iAttVal != attValues->end(); ++iAttVal) {
	        map<G4String,G4AttDef>::const_iterator iAttDef =
	        attDefs->find(iAttVal->GetName());
	        if (iAttDef == attDefs->end()) {
//	            cout << "G4HepRepFileSceneHandler::AddThis(traj):"
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
*/
}

bool G4HepRepSceneHandler::IsEventData () {
    return !currentPV || fReadyForTransients;
}
