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
          geometryLayer         ("Geometry"),
          eventLayer            ("Event"),
          calHitLayer           ("CalHit"),
          trajectoryLayer       ("Trajectory"),
          trajectoryPointLayer  ("TrajectoryPoint"),
          hitLayer              ("Hit"),
          rootVolumeName        ("Geometry"),
          eventNumber           (1),
          currentDepth          (0),
          currentPV             (0),
          currentLV             (0),
          _heprep               (NULL),
          _geometryInstanceTree (NULL),
          _geometryRootInstance (NULL),
          _geometryTypeTree     (NULL),
          _geometryRootType     (NULL),
          _eventInstanceTree    (NULL),
          _eventInstance        (NULL),
          _eventTypeTree        (NULL),
          _eventType            (NULL),
          _trajectoryType       (NULL),
          _trajectoryPointType  (NULL),
          _hitType              (NULL),
          _calHitType           (NULL),
          _calHitFaceType       (NULL)        
{

#ifdef DEBUG
    cout << "G4HepRepSceneHandler::G4HepRepSceneHandler: " << system << endl;
//           << GetScene()->GetName() << endl;
#endif

    materialState[kStateSolid]      = G4String("Solid");
    materialState[kStateLiquid]     = G4String("Liquid");
    materialState[kStateGas]        = G4String("Gas");
    materialState[kStateUndefined]  = G4String("Undefined");    

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

    if (_heprep != NULL) return;

    // all done on demand
}


/**
 * Returns true if the HepRep was (already) closed, false if the HepRep is still open
 */
bool G4HepRepSceneHandler::CloseHepRep() {

#ifdef DEBUG
    cout << "G4HepRepSceneHandler::CloseHepRep() " << endl;
#endif

    if (_heprep == NULL) return true;

    // add geometry to the heprep
    G4HepRepViewer* viewer = dynamic_cast<G4HepRepViewer*>(GetCurrentViewer());
    viewer->ProcessScene();

    // Give this HepRep all of the layer order info for both geometry and event,
    // since these will both end up in a single HepRep.
    if (_geometryRootType    != NULL) _heprep->addLayer(geometryLayer);
    if (_eventType           != NULL) _heprep->addLayer(eventLayer);
    if (_calHitType          != NULL) _heprep->addLayer(calHitLayer);
    if (_trajectoryType      != NULL) _heprep->addLayer(trajectoryLayer);
    if (_trajectoryPointType != NULL) _heprep->addLayer(trajectoryPointLayer);
    if (_hitType             != NULL) _heprep->addLayer(hitLayer);

    // write out the heprep
    char eventName[255];
    sprintf(eventName, "event%010d.heprep", eventNumber);
    writer->write(_heprep, eventName);
    eventNumber++;

    delete _heprep;
    _heprep = NULL;

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

    HepRepInstance* instance = factory->createHepRepInstance(getEventInstance(), getTrajectoryType());

    SetColor(instance, GetColor(line));

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

    HepRepInstance* instance = factory->createHepRepInstance(getEventInstance(), getHitType());

    SetColor(instance, GetColor(line));

    SetMarker(instance, line);

    // Default MarkName is set to Circle for this Type.
    switch (line.GetMarkerType()) {
        case line.dots:
            instance->addAttValue("Fill", true);
            SetColor(instance, GetColor(line), G4String("FillColor"));
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

    HepRepInstance* instance = factory->createHepRepInstance(getEventInstance(), getHitType());

    G4Point3D center = transform * circle.GetPosition();

    SetColor (instance, GetColor(circle));

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

    HepRepInstance* instance;
    if (IsEventData()) {
        instance = factory->createHepRepInstance(getEventInstance(), getCalHitType());
    } else {
        instance = getGeometryInstance(currentLV, currentDepth);
    }
    
    G4bool notLastFace;
    do {
        HepRepInstance* face;
        if (IsEventData()) {
            face = factory->createHepRepInstance(instance, getCalHitType());
        } else {
            face = getGeometryInstance("*Face", currentDepth+1);
        }
        
        SetLine(face, polyhedron);
        SetColor(face, GetColor(polyhedron));
        if (IsEventData()) SetColor(face, GetColor(polyhedron), G4String("FillColor"));

        notLastFace = polyhedron.GetNextNormal (surfaceNormal);

        G4int edgeFlag = 1;
        G4bool notLastEdge;
        do {
	        notLastEdge = polyhedron.GetNextVertex (vertex, edgeFlag);
            vertex = transform * vertex;
            factory->createHepRepPoint(face, vertex.x(), vertex.y(), vertex.z());
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

    HepRepInstance* instance = factory->createHepRepInstance(getEventInstance(), getHitType());

    G4Point3D center = transform * square.GetPosition();

    SetColor (instance, GetColor(square));

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


void G4HepRepSceneHandler::AddThis (const G4VTrajectory& trajectory) {
//#ifdef PDEBUG
    cout << "G4HepRepSceneHandler::AddThis(G4VTrajectory&) " << endl;
//#endif

    vector<G4AttValue>* trajectoryAttValues = trajectory.CreateAttValues();
    const map<G4String,G4AttDef>* trajectoryAttDefs = trajectory.GetAttDefs();

    HepRepInstance* trajectoryInstance = factory->createHepRepInstance(getEventInstanceTree(), getTrajectoryType(trajectoryAttDefs));
    addAttVals(trajectoryInstance, trajectoryAttValues);
    delete trajectoryAttValues;

    // Set the LineColor attribute according to the particle charge.
    float red = 0.;
    float green = 0.;
    float blue = 0.;
    const G4double charge = trajectory.GetCharge();
    if(charge>0.0)      blue =  1.0; // Blue = positive.
    else if(charge<0.0) red  =  1.0; // Red = negative.
    else                green = 1.0; // Green = neutral.
    SetColor(trajectoryInstance, G4Color(red, green, blue));

    // Specify the polyline by using the trajectory points.
    G4int i;
    for (i = 0; i < trajectory.GetPointEntries(); i++) {
        G4VTrajectoryPoint* trajectoryPoint = trajectory.GetPoint(i);
        G4Point3D vertex = trajectoryPoint->GetPosition();
        factory->createHepRepPoint(trajectoryInstance, vertex.x(), vertex.y(), vertex.z());

        // add all points also as separate instances (MD -> JP: why?), with single points
        vector<G4AttValue>* pointAttValues = trajectoryPoint->CreateAttValues();
        const map<G4String,G4AttDef>* pointAttDefs = trajectoryPoint->GetAttDefs();

        HepRepInstance* pointInstance = factory->createHepRepInstance(trajectoryInstance, getTrajectoryPointType(pointAttDefs));
        addAttVals(pointInstance, pointAttValues);
// FIXME seems to crash...
//        delete pointAttValues;

        // Specify the position of the trajectory point.
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
}


void G4HepRepSceneHandler::EndPrimitives () {
#ifdef DEBUG
    cout << "G4HepRepSceneHandler::EndPrimitives" << endl;
#endif
    G4VSceneHandler::EndPrimitives ();
}


void G4HepRepSceneHandler::SetColor (HepRepAttribute *attribute,
				      const G4Color& color,
				      const G4String& key) {
#ifdef VDEBUG
    cout << "G4HepRepSceneHandler::SetColor : red : " << color.GetRed ()   <<
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
        SetColor(instance, GetColor(marker), G4String("FillColor"));
    }
}


/*
HepRepInstance* G4HepRepSceneHandler::CreateGeometryInstance(G4String typeName, G4int depth) {

//    cout << "GeometryInstance: " << typeName << " " << depth << endl;

    return factory->createHepRepInstance(geometryInstanceTree, detectorType);
*/
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
}
*/

void G4HepRepSceneHandler::addAttDefs(HepRepType* type, const map<G4String,G4AttDef>* attDefs) {
    if (attDefs == NULL) return;

    // Specify additional attribute definitions.  Assumes that all attributes for
    //  the given type can be found from the first instance of this type.
    map<G4String,G4AttDef>::const_iterator attDefIterator = attDefs->begin();
    while (attDefIterator != attDefs->end()) {
	    // Protect against incorrect use of Category.  Anything value other than the
	    // standard ones will be considered to be in the physics category.
	    G4String category = attDefIterator->second.GetCategory();
	    if ((category == "Draw") ||
    		(category == "Physics") ||
    		(category == "Association") ||
    		(category == "PickAction")) {

            category = "Physics";
        }
        type->addAttDef(attDefIterator->first, attDefIterator->second.GetDesc(),
                        category, attDefIterator->second.GetExtra());
        attDefIterator++;
    }
}

void G4HepRepSceneHandler::addAttVals(HepRepInstance* instance, vector<G4AttValue>* attValues) {
    if (attValues == NULL) return;

    // Copy the instance's G4AttValues to HepRepAttValues.
    vector<G4AttValue>::iterator attValIterator;
    for (attValIterator = attValues->begin(); attValIterator != attValues->end(); attValIterator++) {
        // Use GetDesc rather than GetName once WIRED can handle names with spaces in them.
        //instance->addAttValue(iAttDef->second.GetDesc(), iAttVal->GetValue());
        instance->addAttValue(attValIterator->GetName(), attValIterator->GetValue());
    }
    delete attValues;  // AttValues must be deleted after use.
}


bool G4HepRepSceneHandler::IsEventData () {
    return !currentPV || fReadyForTransients;
}

HepRep* G4HepRepSceneHandler::getHepRep() {
    if (_heprep == NULL) {
        // Create the HepRep that holds the Trees.
        _heprep = factory->createHepRep();
    }
    return _heprep;
}   

HepRepInstanceTree* G4HepRepSceneHandler::getGeometryInstanceTree() {
    if (_geometryInstanceTree == NULL) {
        // Create the Goemetry InstanceTree.
        _geometryInstanceTree = factory->createHepRepInstanceTree("G4GeometryData", "1.0", getGeometryTypeTree());
        getHepRep()->addInstanceTree(_geometryInstanceTree);
    }
    return _geometryInstanceTree;
}

HepRepInstance* G4HepRepSceneHandler::getGeometryRootInstance() {
    if (_geometryRootInstance == NULL) {
        // Create the top level Geometry Instance.
        _geometryRootInstance = factory->createHepRepInstance(getGeometryInstanceTree(), getGeometryRootType());
    }
    return _geometryRootInstance;
}

HepRepInstance* G4HepRepSceneHandler::getGeometryInstance(G4LogicalVolume* volume, int depth) {
    HepRepInstance* instance = getGeometryInstance(volume->GetName(), depth);

    // FIXME we could avoid some copying here    
    instance->addAttValue("LVol",     volume->GetName());
    instance->addAttValue("Solid",    volume->GetSolid()->GetName());
    instance->addAttValue("EType",    volume->GetSolid()->GetEntityType());
    instance->addAttValue("Material", volume->GetMaterial()->GetName());
    instance->addAttValue("Density",  volume->GetMaterial()->GetDensity());
    instance->addAttValue("Radlen",   volume->GetMaterial()->GetRadlen());
    
    G4String state = materialState[volume->GetMaterial()->GetState()];
    instance->addAttValue("State", state);
    
    return instance;
}

HepRepInstance* G4HepRepSceneHandler::getGeometryInstance(G4String volumeName, int depth) {
    // no extra checks since these are done in the geometryType already

    // adjust depth, also pop the current instance
    while (_geometryInstance.size() > depth) {
        _geometryInstance.pop();
    }
    
    // get parent
    HepRepInstance* parent = (_geometryInstance.empty()) ? getGeometryRootInstance() : _geometryInstance.top();
    
    // get type
    HepRepType* type = getGeometryType(volumeName, depth);
    type->addAttValue("DrawAs", G4String("Polygon"));
    
    // create instance
    HepRepInstance* instance = factory->createHepRepInstance(parent, type);
    _geometryInstance.push(instance);

    return instance;
}

HepRepTypeTree* G4HepRepSceneHandler::getGeometryTypeTree() {
    if (_geometryTypeTree == NULL) {
        // Create the Geometry TypeTree.
        HepRepTreeID* geometryTreeID = factory->createHepRepTreeID("G4GeometryTypes", "1.0");
        _geometryTypeTree = factory->createHepRepTypeTree(geometryTreeID);
        getHepRep()->addTypeTree(_geometryTypeTree);
    }
    return _geometryTypeTree;
}

HepRepType* G4HepRepSceneHandler::getGeometryRootType() {
    if (_geometryRootType == NULL) {
        // Create the top level Geometry Type.
        _geometryRootType = factory->createHepRepType(NULL, rootVolumeName);
        _geometryRootType->addAttValue("Layer", geometryLayer);
    
        // Add attdefs used by all geometry types.
        _geometryRootType->addAttDef("LVol", "Logical Volume", "Physics","");
        _geometryRootType->addAttDef("Solid", "Solid Name", "Physics","");
        _geometryRootType->addAttDef("EType", "Entity Type", "Physics","");
        _geometryRootType->addAttDef("Material", "Material Name", "Physics","");
        _geometryRootType->addAttDef("Density", "Material Density", "Physics","");
        _geometryRootType->addAttDef("State", "Material State", "Physics","");
        _geometryRootType->addAttDef("Radlen", "Material Radiation Length", "Physics","");
        getGeometryTypeTree()->addType(_geometryRootType);
        
        _geometryType["/"+_geometryRootType->getName()] = _geometryRootType;
    }
    return _geometryRootType;
}

HepRepType* G4HepRepSceneHandler::getGeometryType(G4String volumeName, int depth) {
    // make sure we have a root
    getGeometryRootType();   

    // construct the full name for this volume
    G4String name = getFullTypeName(volumeName, depth);
    
    // lookup type and create if necessary
    HepRepType* type = _geometryType[name];    
    if (type == NULL) {
        G4String parentName = getParentTypeName(depth);
        HepRepType* parentType = _geometryType[parentName];
        // FIXME should have short path element name...(volumeName)
        type = factory->createHepRepType(parentType, name);
        _geometryType[name] = type;
    }
    return type;   
}

G4String G4HepRepSceneHandler::getFullTypeName(G4String volumeName, int depth) {
    // check for name depth
    if (depth > _geometryTypeName.size()) {
        // there is a problem, book this type under problems
        G4String problem = "HierarchyProblem";
        if (_geometryType["/"+problem] == NULL) {
            // FIXME should have short path element name...
            HepRepType* type = factory->createHepRepType(getGeometryRootType(), "/"+problem);
            _geometryType["/"+problem] = type;
        }
        return "/" + problem + "/" + volumeName;
    }
    
    // adjust name depth, also pop the current volumeName
    while (_geometryTypeName.size() > depth) {
        _geometryTypeName.pop_back();
    }
    
    // construct full name and push it
    G4String name = (_geometryTypeName.empty()) ? "/"+rootVolumeName : _geometryTypeName.back();
    name = name + "/" + volumeName;
    _geometryTypeName.push_back(name);
    return name;
}

G4String G4HepRepSceneHandler::getParentTypeName(int depth) {
    return (depth >= 1) ? _geometryTypeName[depth-1] : "/"+rootVolumeName;
}       

HepRepInstanceTree* G4HepRepSceneHandler::getEventInstanceTree() {
    if (_eventInstanceTree == NULL) {
        // Create the Event InstanceTree.
        _eventInstanceTree = factory->createHepRepInstanceTree("G4EventData", "1.0", getEventTypeTree());
        getHepRep()->addInstanceTree(_eventInstanceTree);
    }
    return _eventInstanceTree;
}

HepRepInstance* G4HepRepSceneHandler::getEventInstance() {
    if (_eventInstance == NULL) {
        // Create the top level Event Instance.
        _eventInstance = factory->createHepRepInstance(getEventInstanceTree(), getEventType());
    }
    return _eventInstance;
}

HepRepTypeTree* G4HepRepSceneHandler::getEventTypeTree() {
    if (_eventTypeTree == NULL) {
        // Create the Event TypeTree.
        HepRepTreeID* eventTreeID = factory->createHepRepTreeID("G4EventTypes", "1.0");
        _eventTypeTree = factory->createHepRepTypeTree(eventTreeID);
        getHepRep()->addTypeTree(_eventTypeTree);
    }
    return _eventTypeTree;
}

HepRepType* G4HepRepSceneHandler::getEventType() {
    if (_eventType == NULL) {
        // Create the top level Event Type.
        _eventType = factory->createHepRepType(NULL, "Event");
        _eventType->addAttValue("Layer", eventLayer);
    
        getEventTypeTree()->addType(_eventType);
    }
    return _eventType;
}

HepRepType* G4HepRepSceneHandler::getTrajectoryType(const map<G4String,G4AttDef>* attDefs) {
    if (_trajectoryType == NULL) {
        _trajectoryType = factory->createHepRepType(getEventType(), "Trajectory");
        _trajectoryType->addAttValue("Layer", trajectoryLayer);
        _trajectoryType->addAttValue("DrawAs", "Line");

        addAttDefs(_trajectoryType, attDefs);
    }
    return _trajectoryType;
}

HepRepType* G4HepRepSceneHandler::getTrajectoryPointType(const map<G4String,G4AttDef>* attDefs) {
    if (_trajectoryPointType == NULL) {
        _trajectoryPointType = factory->createHepRepType(getTrajectoryType(), "Trajectory/Point");
        _trajectoryPointType->addAttValue("Layer", trajectoryPointLayer);
        _trajectoryPointType->addAttValue("DrawAs", "Circle");
        _trajectoryPointType->addAttValue("MarkName", "Square");
        _trajectoryPointType->addAttValue("MarkType", "Real");
        _trajectoryPointType->addAttValue("Fill", true);
        _trajectoryPointType->addAttValue("Visibility", true);

        addAttDefs(_trajectoryPointType, attDefs);
    }
    return _trajectoryPointType;
}

HepRepType* G4HepRepSceneHandler::getHitType(const map<G4String,G4AttDef>* attDefs) {
    if (_hitType == NULL) {
        _hitType = factory->createHepRepType(getEventType(), "Hit");
        _hitType->addAttValue("Layer", hitLayer);
        _hitType->addAttValue("DrawAs", "Point");
        _hitType->addAttValue("MarkName", "Circle");
        _hitType->addAttValue("MarkType", "Real");
        _hitType->addAttValue("Fill", true);

        addAttDefs(_hitType, attDefs);
    }
    return _hitType;
}

HepRepType* G4HepRepSceneHandler::getCalHitType(const map<G4String,G4AttDef>* attDefs) {
    if (_calHitType == NULL) {
        _calHitType = factory->createHepRepType(getEventType(), "CalHit");
        _calHitType->addAttValue("Layer", calHitLayer);
        _calHitType->addAttValue("Fill", true);

        addAttDefs(_calHitType, attDefs);
    }
    return _calHitType;
}

HepRepType* G4HepRepSceneHandler::getCalHitFaceType(const map<G4String,G4AttDef>* attDefs) {
    if (_calHitFaceType == NULL) {
        _calHitFaceType = factory->createHepRepType(getCalHitType(), "CalHit/Face");
        _calHitFaceType->addAttValue("Layer", calHitLayer);
        _calHitFaceType->addAttValue("DrawAs", "Line");

        addAttDefs(_calHitFaceType, attDefs);
    }
    return _calHitFaceType;
}

