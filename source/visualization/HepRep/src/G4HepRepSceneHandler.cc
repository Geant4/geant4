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
// NOTE not available on Solaris 5.2 and Linux g++ 2.95.2
// #include <sstream>
#include <iomanip>
#include <fstream>
#include <cmath>

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

//#define SDEBUG 1
//#define PDEBUG 1

G4HepRepSceneHandler::G4HepRepSceneHandler (G4VGraphicsSystem& system, G4HepRepMessenger& heprepMessenger, const G4String& name)
        : G4VSceneHandler (system, sceneCount++, name),
          messenger             (heprepMessenger),
          geometryLayer         ("Geometry"),
          eventLayer            ("Event"),
          calHitLayer           ("CalHit"),
          trajectoryLayer       ("Trajectory"),
          hitLayer              ("Hit"),
          rootVolumeName        ("Geometry"),
          baseName              (""),
          eventNumberPrefix     (""),
          eventNumberSuffix     (""),
          eventNumber           (1),
          eventNumberWidth      (-1),
          extension             (""),
          writeMultipleFiles    (false),
          currentDepth          (0),
          currentPV             (0),
          currentLV             (0),
          _heprep               (NULL)
{

#ifdef SDEBUG
    cout << "G4HepRepSceneHandler::G4HepRepSceneHandler: " << system << endl;
#endif

    materialState[kStateSolid]      = G4String("Solid");
    materialState[kStateLiquid]     = G4String("Liquid");
    materialState[kStateGas]        = G4String("Gas");
    materialState[kStateUndefined]  = G4String("Undefined");    

    factory = new XMLHepRepFactory();
    writer = NULL;
    
    // opening of file deferred to closeHepRep();
    openHepRep();
}


G4HepRepSceneHandler::~G4HepRepSceneHandler () {
#ifdef SDEBUG
    cout << "G4HepRepSceneHandler::~G4HepRepSceneHandler() " << endl;
#endif
    close();

    delete factory;
    factory = NULL;

    dynamic_cast<G4HepRep*>(GetGraphicsSystem())->removeSceneHandler();
}


void G4HepRepSceneHandler::open(G4String name) {
    if (writer != NULL) return;

    if (name == "stdout") {
#ifdef SDEBUG
        cout << "G4HepRepSceneHandler::Open() stdout" << endl;
#endif
        writer = factory->createHepRepWriter(&cout, false, false);
        out  = NULL;
        baseName = name;
        eventNumberPrefix = "";
        eventNumberSuffix = "";
        extension = "";
        writeMultipleFiles = false;
        eventNumber = 0;
        eventNumberWidth = 0;        
    } else if (name == "stderr") {
#ifdef SDEBUG
        cout << "G4HepRepSceneHandler::Open() stderr" << endl;
#endif
        writer = factory->createHepRepWriter(&cerr, false, false);
        out = NULL;
        baseName = name;
        eventNumberPrefix = "";
        eventNumberSuffix = "";
        extension = "";
        writeMultipleFiles = false;
        eventNumber = 0;
        eventNumberWidth = 0;        
    } else {
#ifdef SDEBUG
        cout << "G4HepRepSceneHandler::Open() " << name << endl;
#endif
        if (eventNumberWidth < 0) {
            // derive filename(s)
            // check for extensions
            const unsigned int numberOfExtensions = 4;
            string ext[numberOfExtensions] = {".heprep", ".heprep.xml", ".heprep.zip", ".heprep.gz"};
            unsigned int i=0;
            while (i < numberOfExtensions) {
                int dot = name.size() - ext[i].size();
                if ((dot >= 0) && 
                    (name.substr(dot, ext[i].size()) == ext[i])) break;
                i++;
            }
        
            if (i != numberOfExtensions) {
                extension = ext[i];
                int dot = name.length() - extension.length();
                baseName = (dot >= 0) ? name.substr(0, dot) : "";
            } else {  
                extension = ".heprep.zip";
                baseName = name;
            }
        
            writeMultipleFiles = false;
            int startDigit = -1; int endDigit = -1;
            string suffix = messenger.getEventNumberSuffix();
            if (suffix != "") {
                // look for 0000 pattern in suffix
                endDigit = suffix.length()-1; 
                while (endDigit >= 0) {
                    if (isdigit(suffix.at(endDigit))) break;
                    endDigit--;
                }
                if (endDigit < 0) {
                    cerr << "/vis/heprep/appendEventNumberSuffix contains no digits" << endl;
                } else {
                    writeMultipleFiles = true;
                    startDigit = endDigit;
                    while (startDigit >= 0) {
                        if (!isdigit(suffix.at(startDigit))) break;
                        startDigit--;
                    }
                    startDigit++;
                }                
            }

            if (writeMultipleFiles) {
                eventNumberPrefix = suffix.substr(0, startDigit);
                eventNumber = atoi(suffix.substr(startDigit, endDigit).c_str());
                eventNumberWidth = endDigit +1 - startDigit;
                eventNumberSuffix = suffix.substr(endDigit+1);
            } else {
                // open single file here
                openFile(baseName+extension);
        
                eventNumber = 1;
                eventNumberWidth = 10;
                eventNumberPrefix = "";
                eventNumberSuffix = "";
            }   
        }    
    }
}


void G4HepRepSceneHandler::openHepRep() {
#ifdef SDEBUG
    cout << "G4HepRepSceneHandler::OpenHepRep() " << endl;
#endif

    if (_heprep != NULL) return;

    // all done on demand, once pointers are set to NULL
    _geometryInstanceTree = NULL;
    _geometryRootInstance = NULL;
    _geometryInstance.clear();
    _geometryTypeTree     = NULL;
    _geometryRootType     = NULL;
    _geometryTypeName.clear();
    _geometryType.clear();
    _eventInstanceTree    = NULL;
    _eventInstance        = NULL;
    _eventTypeTree        = NULL;
    _eventType            = NULL;
    _trajectoryType       = NULL;
    _hitType              = NULL;
    _calHitType           = NULL;
    _calHitFaceType       = NULL;     
}


/**
 * Returns true if the HepRep was (already) closed, false if the HepRep is still open
 */
bool G4HepRepSceneHandler::closeHepRep(bool final) {
    if (_heprep == NULL) return true;

#ifdef SDEBUG
    cout << "G4HepRepSceneHandler::CloseHepRep() start" << endl;
#endif

    // if this is the final close, then there should not be any event pending to be written.
    if (final) {
        if (_eventInstanceTree != NULL) {
            cerr << "WARNING: you probably used '/vis/viewer/endOfEventAction accumulate' and "
                 << "forgot to call /vis/viewer/update before exit. No event written." << endl;
        }
    } else {

        // add geometry to the heprep if there is an event (separate geometries are written
        // using DrawView() called from /vis/viewer/flush)
        if (_eventInstanceTree != NULL) {
            GetCurrentViewer()->DrawView();
    
            // couple geometry to event if geometry was written
            if ((_geometryInstanceTree != NULL)) {
                getEventInstanceTree()->addInstanceTree(getGeometryInstanceTree());
            }
        }
    
        // force inclusion of all subtypes of event
        if (_eventInstanceTree != NULL) {
            getEventType();
            getTrajectoryType();
            getHitType();
            getCalHitType();
            getCalHitFaceType();
        }
    
        // Give this HepRep all of the layer order info for both geometry and event,
        // since these will both end up in a single HepRep.
        if (_geometryRootType    != NULL) _heprep->addLayer(geometryLayer);
        if (_eventType           != NULL) _heprep->addLayer(eventLayer);
        if (_calHitType          != NULL) _heprep->addLayer(calHitLayer);
        if (_trajectoryType      != NULL) _heprep->addLayer(trajectoryLayer);
        if (_hitType             != NULL) _heprep->addLayer(hitLayer);

        // open heprep file
        if (writer == NULL) {
            open((GetScene() == NULL) ? G4String("G4HepRepOutput.heprep.zip") : GetScene()->GetName());
        }
    
        if (writeMultipleFiles) {
// NOTE: does not work on Solaris 5.2 and Linux 2.95.2
//            stringstream fileName;
//            fileName << baseName << eventNumberPrefix << setw(eventNumberWidth) << setfill('0') << eventNumber << eventNumberSuffix << extension;
//            openFile(fileName.str());
// Use instead:
            char fileName[128];
            char fileFormat[128];
            sprintf(fileFormat, "%s%d%s", "%s%s%0", eventNumberWidth, "d%s%s");
            sprintf(fileName, fileFormat, baseName.c_str(), eventNumberPrefix.c_str(), eventNumber, eventNumberSuffix.c_str(), extension.c_str());
            openFile(G4String(fileName));
        }
            
        // write out the heprep
// NOTE: does not work on Solaris 5.2 and Linux 2.95.2
//        stringstream eventName;
//        eventName << "event-" << setw(eventNumberWidth) << setfill('0') << eventNumber << ".heprep";
//        writer->write(_heprep, eventName.str());
// Use instead:
        char eventName[128];
        char eventFormat[128];
        sprintf(eventFormat, "%s%d%s", "event-%0", eventNumberWidth, "d.heprep");
        sprintf(eventName, eventFormat, eventNumber);
        writer->write(_heprep, G4String(eventName));

        eventNumber++;
    }
    
    delete _heprep;
    _heprep = NULL;

    if (writeMultipleFiles) closeFile();

    return true;
}


void G4HepRepSceneHandler::close() {

#ifdef SDEBUG
    cout << "G4HepRepSceneHandler::Close() " << endl;
#endif

    if (writer == NULL) return;

    if (!writeMultipleFiles) {
        closeHepRep(true);
        closeFile();
    }
}

void G4HepRepSceneHandler::openFile(G4String name) {
    out = new ofstream(name.c_str(), std::ios::out | std::ios::binary );
    writer = factory->createHepRepWriter(out, extension == ".heprep.zip", (extension == ".heprep.zip") || (extension == ".heprep.gz"));
}

void G4HepRepSceneHandler::closeFile() {
    writer->close();
    delete writer;
    writer = NULL;

    delete out;
    out = NULL;
}


void G4HepRepSceneHandler::EstablishSpecials(G4PhysicalVolumeModel& model) {
#ifdef SDEBUG
    cout << "G4HepRepSceneHandler::EstablishSpecials(G4PhysicalVolumeModel&) " << endl;
#endif
    model.DefinePointersToWorkingSpace(&currentDepth, &currentPV, &currentLV);
}


void G4HepRepSceneHandler::BeginModeling() {
#ifdef SDEBUG
    cout << "G4HepRepSceneHandler::BeginModeling() " << endl;
#endif
    G4VSceneHandler::BeginModeling();
}


void G4HepRepSceneHandler::EndModeling() {
#ifdef SDEBUG
    cout << "G4HepRepSceneHandler::EndModeling() " << endl;
#endif
    G4VSceneHandler::EndModeling();
}


void G4HepRepSceneHandler::AddPrimitive (const G4Polyline& line) {

#ifdef PDEBUG
    cout << "G4HepRepSceneHandler::AddPrimitive(G4Polyline&) " << line.size() << endl;
#endif

    HepRepInstance* instance = factory->createHepRepInstance(getEventInstance(), getTrajectoryType());

    setColor(instance, GetColor(line));

    setVisibility(instance, line);

    setLine(instance, line);

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

    setColor(instance, GetColor(line));

    setVisibility(instance, line);

    setMarker(instance, line);

    // Default MarkName is set to Circle for this Type.
    switch (line.GetMarkerType()) {
        case line.dots:
            setAttribute(instance, "Fill", true);
            setColor(instance, GetColor(line), G4String("FillColor"));
            break;
        case line.circles:
            break;
        case line.squares:
            setAttribute(instance, "MarkName", G4String("Box"));
            break;
        case line.line:
        default:
            setAttribute(instance, "MarkName", G4String("Plus"));
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

    setColor (instance, GetColor(circle));

    setVisibility(instance, circle);

    setMarker(instance, circle);

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
    if (isEventData()) {
        instance = factory->createHepRepInstance(getEventInstance(), getCalHitType());
    } else {
        instance = getGeometryInstance(currentLV, currentDepth);
    }
        
    setVisibility(instance, polyhedron);
	
    G4bool notLastFace;
    do {
        HepRepInstance* face;
        if (isEventData()) {
            face = factory->createHepRepInstance(instance, getCalHitFaceType());
        } else {
            face = getGeometryInstance("*Face", currentDepth+1);
            setAttribute(face, "PickParent", true);
        }
        
        setLine(face, polyhedron);
        setColor(face, GetColor(polyhedron));
        if (isEventData()) setColor(face, GetColor(polyhedron), G4String("FillColor"));

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

    setColor (instance, GetColor(square));

    setVisibility(instance, square);

    setMarker(instance, square);

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
#ifdef PDEBUG
    cout << "G4HepRepSceneHandler::AddThis(G4VTrajectory&) " << endl;
#endif

    vector<G4AttValue>* trajectoryAttValues = trajectory.CreateAttValues();
    const map<G4String,G4AttDef>* trajectoryAttDefs = trajectory.GetAttDefs();

    HepRepType* trajectoryType = getTrajectoryType();
    addAttDefs(trajectoryType, trajectoryAttDefs);
    
    // these attValues are non-standard, so can only be added when we have the attDef.
    double charge = 0;
    _trajectoryType->addAttValue("Ch", charge);
    G4Color color = getColor(charge);
    _trajectoryType->addAttValue("Color", color.GetRed(), color.GetGreen(), color.GetBlue(), color.GetAlpha());
    _trajectoryType->addAttValue("ID", -1);
    _trajectoryType->addAttValue("IMom", G4String(""));
    _trajectoryType->addAttValue("PDG", -1);
    _trajectoryType->addAttValue("PN", G4String(""));
    _trajectoryType->addAttValue("PID", -1);

    HepRepInstance* trajectoryInstance = factory->createHepRepInstance(getEventInstance(), trajectoryType);
    addAttVals(trajectoryInstance, trajectoryAttDefs, trajectoryAttValues);
    
    delete trajectoryAttValues;

    setColor(trajectoryInstance, getColor(trajectory.GetCharge()));    
    setAttribute(trajectoryInstance, "LineWidth", 1.0);

    // Specify the polyline by using the trajectory points.
    for (int i = 0; i < trajectory.GetPointEntries(); i++) {
        G4VTrajectoryPoint* trajectoryPoint = trajectory.GetPoint(i);
        G4Point3D vertex = trajectoryPoint->GetPosition();
        HepRepPoint* point = factory->createHepRepPoint(trajectoryInstance, vertex.x(), vertex.y(), vertex.z());


        if (messenger.addPointAttributes()) {
            vector<G4AttValue>* pointAttValues = trajectoryPoint->CreateAttValues();
            const map<G4String,G4AttDef>* pointAttDefs = trajectoryPoint->GetAttDefs();
            addAttVals(point, pointAttDefs, pointAttValues);
            delete pointAttValues;
        }
    }
}



void G4HepRepSceneHandler::PreAddThis (const G4Transform3D& objectTransformation,
				  const G4VisAttributes& visAttribs) {

    G4VSceneHandler::PreAddThis (objectTransformation, visAttribs);

    transform = objectTransformation;
#ifdef SDEBUG
    cout << "G4HepRepSceneHandler::PreAddThis(G4Transform3D&, G4VisAttributes&)" << endl;
#endif
}


void G4HepRepSceneHandler::PostAddThis () {
#ifdef SDEBUG
    cout << "G4HepRepSceneHandler::PostAddThis()" << endl;
#endif
    G4VSceneHandler::PostAddThis();
}


void G4HepRepSceneHandler::BeginPrimitives (const G4Transform3D& objectTransformation) {
#ifdef SDEBUG
    cout << "G4HepRepSceneHandler::BeginPrimitives(G4Transform3D&)" << endl;
#endif

    G4VSceneHandler::BeginPrimitives (objectTransformation);
    transform = objectTransformation;
}


void G4HepRepSceneHandler::EndPrimitives () {
#ifdef SDEBUG
    cout << "G4HepRepSceneHandler::EndPrimitives" << endl;
#endif
    G4VSceneHandler::EndPrimitives ();
}


void G4HepRepSceneHandler::setColor (HepRepAttribute *attribute,
				      const G4Color& color,
				      const G4String& key) {
#ifdef CDEBUG
    cout << "G4HepRepSceneHandler::setColor : red : " << color.GetRed ()   <<
                                  " green : " << color.GetGreen () <<
                                  " blue : " << color.GetBlue ()   << endl;
#endif

    setAttribute(attribute, key, color.GetRed(), color.GetGreen(), color.GetBlue(),color.GetAlpha());
}

G4Color G4HepRepSceneHandler::getColor (G4double charge) {
    float red = 0.;
    float green = 0.;
    float blue = 0.;
    if(charge>0.0)      blue =  1.0; // Blue = positive.
    else if(charge<0.0) red  =  1.0; // Red = negative.
    else                green = 1.0; // Green = neutral.
    return G4Color(red, green, blue);
}

void G4HepRepSceneHandler::setVisibility (HepRepAttribute *attribute, const G4Visible& visible) {
    const G4VisAttributes* atts = visible.GetVisAttributes();

    setAttribute(attribute, "Visibility", (atts && (atts->IsVisible()==0)) ? false : true);
}

void G4HepRepSceneHandler::setLine (HepRepAttribute *attribute, const G4Visible& visible) {
    const G4VisAttributes* atts = visible.GetVisAttributes();
    
    setAttribute(attribute, "LineWidth", (atts != NULL) ? atts->GetLineWidth() : 1.0);
    
    if (atts != NULL) {
        switch (atts->GetLineStyle()) {
            case G4VisAttributes::dotted:
                setAttribute(attribute, "LineStyle", G4String("Dotted"));
                break;
            case G4VisAttributes::dashed:
                setAttribute(attribute, "LineStyle", G4String("Dashed"));
                break;
            case G4VisAttributes::unbroken:
            default:
                break;
        }
    }
}

void G4HepRepSceneHandler::setMarker (HepRepAttribute *attribute, const G4VMarker& marker) {
    MarkerSizeType markerType;
    G4double size = GetMarkerRadius( marker , markerType );

    setAttribute(attribute, "MarkSize", size);

    if (markerType == screen) setAttribute(attribute, "MarkType", G4String("Symbol"));
    if (marker.GetFillStyle() == G4VMarker::noFill) {
        setAttribute(attribute, "Fill", false);
    } else {
        setColor(attribute, GetColor(marker), G4String("FillColor"));
    }
}

void G4HepRepSceneHandler::setAttribute(HepRepAttribute* attribute, G4String name, G4String value) {
    HepRepAttValue* attValue = attribute->getAttValue(name);
    if ((attValue == NULL) || (attValue->getString() != value)) {
        HepRepPoint* point = dynamic_cast<HepRepPoint*>(attribute);
        if (point != NULL) {
            if (point->getInstance()->getAttValueFromNode(name) == NULL) {
                attribute = point->getInstance();
            }
        }
        
        HepRepInstance* instance = dynamic_cast<HepRepInstance*>(attribute);
        if (instance != NULL) {
            // look for definition on type (node only)
            if (instance->getType()->getAttValueFromNode(name) == NULL) {
                attribute = instance->getType();
            }   
        }
        
        attribute->addAttValue(name, value);
    }
}

void G4HepRepSceneHandler::setAttribute(HepRepAttribute* attribute, G4String name, bool value) {
    HepRepAttValue* attValue = attribute->getAttValue(name);
    if ((attValue == NULL) || (attValue->getBoolean() != value)) {
        HepRepPoint* point = dynamic_cast<HepRepPoint*>(attribute);
        if (point != NULL) {
            if (point->getInstance()->getAttValueFromNode(name) == NULL) {
                attribute = point->getInstance();
            }
        }
        
        HepRepInstance* instance = dynamic_cast<HepRepInstance*>(attribute);
        if (instance != NULL) {
            // look for definition on type (node only)
            if (instance->getType()->getAttValueFromNode(name) == NULL) {
                attribute = instance->getType();
            }   
        }
        
        attribute->addAttValue(name, value);
    }
}

void G4HepRepSceneHandler::setAttribute(HepRepAttribute* attribute, G4String name, double value) {
    HepRepAttValue* attValue = attribute->getAttValue(name);
    if ((attValue == NULL) || (attValue->getDouble() != value)) {
        HepRepPoint* point = dynamic_cast<HepRepPoint*>(attribute);
        if (point != NULL) {
            if (point->getInstance()->getAttValueFromNode(name) == NULL) {
                attribute = point->getInstance();
            }
        }
        
        HepRepInstance* instance = dynamic_cast<HepRepInstance*>(attribute);
        if (instance != NULL) {
            // look for definition on type (node only)
            if (instance->getType()->getAttValueFromNode(name) == NULL) {
                attribute = instance->getType();
            }   
        }
        
        attribute->addAttValue(name, value);
    }
}

void G4HepRepSceneHandler::setAttribute(HepRepAttribute* attribute, G4String name, int value) {
    HepRepAttValue* attValue = attribute->getAttValue(name);
    if ((attValue == NULL) || (attValue->getInteger() != value)) {
        HepRepPoint* point = dynamic_cast<HepRepPoint*>(attribute);
        if (point != NULL) {
            if (point->getInstance()->getAttValueFromNode(name) == NULL) {
                attribute = point->getInstance();
            }
        }
        
        HepRepInstance* instance = dynamic_cast<HepRepInstance*>(attribute);
        if (instance != NULL) {
            // look for definition on type (node only)
            if (instance->getType()->getAttValueFromNode(name) == NULL) {
                attribute = instance->getType();
            }   
        }
        
        attribute->addAttValue(name, value);
    }
}

void G4HepRepSceneHandler::setAttribute(HepRepAttribute* attribute, G4String name, double red, double green, double blue, double alpha) {
    HepRepAttValue* attValue = attribute->getAttValue(name);
    vector<double> color;
    if (attValue != NULL) color = attValue->getColor();
    if ((color.size() == 0) ||
        (color[0] != red) ||
        (color[1] != green) ||
        (color[2] != blue) ||
        ((color.size() > 3) && (color[3] != alpha))) {
        
        HepRepPoint* point = dynamic_cast<HepRepPoint*>(attribute);
        if (point != NULL) {
            if (point->getInstance()->getAttValueFromNode(name) == NULL) {
                attribute = point->getInstance();
            }
        }
        
        HepRepInstance* instance = dynamic_cast<HepRepInstance*>(attribute);
        if (instance != NULL) {
            // look for definition on type (node only)
            if (instance->getType()->getAttValueFromNode(name) == NULL) {
                attribute = instance->getType();
            }   
        }
        
        attribute->addAttValue(name, red, green, blue, alpha);
    }
}

void G4HepRepSceneHandler::addAttDefs(HepRepDefinition* definition, const map<G4String,G4AttDef>* attDefs) {
    if (attDefs == NULL) return;

    // Specify additional attribute definitions.
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
        definition->addAttDef(attDefIterator->first, attDefIterator->second.GetDesc(),
                        category, attDefIterator->second.GetExtra());
        attDefIterator++;
    }
}

void G4HepRepSceneHandler::addAttVals(HepRepAttribute* attribute, const map<G4String,G4AttDef>* attDefs, vector<G4AttValue>* attValues) {
    if (attValues == NULL) return;

    // Copy the instance's G4AttValues to HepRepAttValues.
    for (vector<G4AttValue>::iterator attValIterator = attValues->begin(); attValIterator != attValues->end(); attValIterator++) {
        // Use GetDesc rather than GetName once WIRED can handle names with spaces in them.
        //attribute->addAttValue(iAttDef->second.GetDesc(), iAttVal->GetValue());
        
        G4String name = attValIterator->GetName();

        HepRepPoint* point = dynamic_cast<HepRepPoint*>(attribute);
        if ((name == "Pos") && (point != NULL)) {
            G4String pos = attValIterator->GetValue();
//            cout << "Pos* " << pos << endl;
            int s = 0;
            int n = 0;
            int m = 0;
            G4String unit;
            for (unsigned int i=0; i<pos.length(); i++) {
                if (pos[i] == ' ') {
                    if (n == 0) {
                        // first coordinate
                        double factor = atof(pos.substr(s, i-s).c_str())/point->getX();
                        m = (int)(log10(factor)+((factor < 1) ? -0.5 : 0.5));
//                        cout << factor << ", " << m << endl;
                    } else if (n == 3) {
                        // unit
                        unit = pos.substr(s, i-s);
                        if (unit == G4String("mum")) {
                            m += -6;
                        } else if (unit == G4String("mm")) {
                            m += -3;
                        } else if (unit == G4String("cm")) {
                            m += -2;
                        } else if (unit == G4String("m")) {
                            m += 0;
                        } else if (unit == G4String("km")) {
                            m += 3;
                        } else {
                            cerr << "HepRepSceneHandler: Unrecognized Unit: '" << unit << "'" << endl;
                        }
                    }
                    s = i+1;
                    n++;
                }
            }
            switch(m) {
                case -6:
                    unit = G4String("mum");
                    break;
                case -3:
                    unit = G4String("mm");
                    break;
                case -2:
                    unit = G4String("cm");
                    break;
                case 0:
                    unit = G4String("m");
                    break;
                case 3:
                    unit = G4String("km");
                    break;
                default:
                    cerr << "HepRepSceneHandler: No valid unit found for m: " << m << endl;
                    unit = G4String("*m");                        
                    break;
            }
//            cout << "U: " << unit << endl;
            setAttribute(attribute, G4String("PointUnit"), unit);
            continue;
        }
                
        // NTP already in points being written
        if (name == "NTP") continue;
        
        // find type of attribute using def
        const map<G4String,G4AttDef>::const_iterator attDefIterator = attDefs->find(name);
        G4String type = attDefIterator->second.GetValueType();
                
        // set based on type
        if ((type == "G4double") || (type == "double")) {   
            setAttribute(attribute, attValIterator->GetName(), atof(attValIterator->GetValue()));           
        } else if ((type == "G4int") || (type == "int")) {  
            setAttribute(attribute, attValIterator->GetName(), atoi(attValIterator->GetValue()));           
        } else { // G4String, string and others
            setAttribute(attribute, attValIterator->GetName(), attValIterator->GetValue());
        }
    }
}


bool G4HepRepSceneHandler::isEventData () {
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
        // Create the Geometry InstanceTree.
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

    setAttribute(instance, "LVol",     volume->GetName());
    setAttribute(instance, "Solid",    volume->GetSolid()->GetName());
    setAttribute(instance, "EType",    volume->GetSolid()->GetEntityType());
    setAttribute(instance, "Material", volume->GetMaterial()->GetName());
    setAttribute(instance, "Density",  volume->GetMaterial()->GetDensity());
    setAttribute(instance, "Radlen",   volume->GetMaterial()->GetRadlen());
    
    G4String state = materialState[volume->GetMaterial()->GetState()];
    setAttribute(instance, "State", state);
    
    return instance;
}

HepRepInstance* G4HepRepSceneHandler::getGeometryInstance(G4String volumeName, int depth) {
    // no extra checks since these are done in the geometryType already

    // adjust depth, also pop the current instance
    while ((int)_geometryInstance.size() > depth) {
        _geometryInstance.pop_back();
    }
    
    // get parent
    HepRepInstance* parent = (_geometryInstance.empty()) ? getGeometryRootInstance() : _geometryInstance.back();
    
    // get type
    HepRepType* type = getGeometryType(volumeName, depth);
    
    // create instance
    HepRepInstance* instance = factory->createHepRepInstance(parent, type);
    _geometryInstance.push_back(instance);

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
        _geometryRootType->addAttDef  ("LVol", "Logical Volume", "Physics","");
        _geometryRootType->addAttValue("LVol", G4String(""));
        _geometryRootType->addAttDef  ("Solid", "Solid Name", "Physics","");
        _geometryRootType->addAttValue("Solid", G4String(""));
        _geometryRootType->addAttDef  ("EType", "Entity Type", "Physics","");
        _geometryRootType->addAttValue("EType", G4String("G4Box"));
        _geometryRootType->addAttDef  ("Material", "Material Name", "Physics","");
        _geometryRootType->addAttValue("Material", G4String("Air"));
        _geometryRootType->addAttDef  ("Density", "Material Density", "Physics","");
        _geometryRootType->addAttValue("Density", 0.0);
        _geometryRootType->addAttDef  ("State", "Material State", "Physics","");
        _geometryRootType->addAttValue("State", G4String("Gas"));
        _geometryRootType->addAttDef  ("Radlen", "Material Radiation Length", "Physics","");
        _geometryRootType->addAttValue("Radlen", 0.0);
        
        // add defaults for Geometry
        _geometryRootType->addAttValue("Color", 0.8, 0.8, 0.8, 1.0);
        _geometryRootType->addAttValue("Visibility", true);
        _geometryRootType->addAttValue("FillColor", 0.8, 0.8, 0.8, 1.0);
        _geometryRootType->addAttValue("LineWidth", 1.0);
        _geometryRootType->addAttValue("DrawAs", G4String("Polygon"));
        _geometryRootType->addAttValue("PickParent", false);
        _geometryRootType->addAttValue("ShowParentAttributes", true);
        
        _geometryRootType->addAttValue("MarkSizeMultiplier", 4.0);
        _geometryRootType->addAttValue("LineWidthMultiplier", 1.0);
        
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
        // HepRep uses hierarchical names
        type = factory->createHepRepType(parentType, volumeName);
        _geometryType[name] = type;
    }
    return type;   
}

G4String G4HepRepSceneHandler::getFullTypeName(G4String volumeName, int depth) {
    // check for name depth
    if (depth > (int)_geometryTypeName.size()) {
        // there is a problem, book this type under problems
        G4String problem = "HierarchyProblem";
        if (_geometryType["/"+problem] == NULL) {
            // HepRep uses hierarchical names
            HepRepType* type = factory->createHepRepType(getGeometryRootType(), problem);
            _geometryType["/"+problem] = type;
        }
        return "/" + problem + "/" + volumeName;
    }
    
    // adjust name depth, also pop the current volumeName
    while ((int)_geometryTypeName.size() > depth) {
        _geometryTypeName.pop_back();
    }
    
    // construct full name and push it
    G4String name = (_geometryTypeName.empty()) ? G4String("/"+rootVolumeName) : _geometryTypeName.back();
    name = name + "/" + volumeName;
    _geometryTypeName.push_back(name);
    return name;
}

G4String G4HepRepSceneHandler::getParentTypeName(int depth) {
    return (depth >= 1) ? _geometryTypeName[depth-1] : G4String("/"+rootVolumeName);
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

        // add defaults for Events
        _eventType->addAttValue("Visibility", true);
        _eventType->addAttValue("Color", 1.0, 1.0, 1.0, 1.0);
        _eventType->addAttValue("FillColor", 1.0, 1.0, 1.0, 1.0);
        _eventType->addAttValue("LineWidth", 1.0);
        _eventType->addAttValue("HasFrame", true);
        _eventType->addAttValue("PickParent", false);
        _eventType->addAttValue("ShowParentAttributes", false);

        _eventType->addAttValue("MarkSizeMultiplier", 4.0);
        _eventType->addAttValue("LineWidthMultiplier", 1.0);

        // Some non-standard attributes
        _eventType->addAttDef("PointUnit", "Length", "Physics", "");
        _eventType->addAttValue("PointUnit", G4String("m"));
            
        getEventTypeTree()->addType(_eventType);
    }
    
    return _eventType;
}

HepRepType* G4HepRepSceneHandler::getTrajectoryType() {
    if (_trajectoryType == NULL) {
        _trajectoryType = factory->createHepRepType(getEventType(), "Trajectory");
        
        _trajectoryType->addAttValue("Layer", trajectoryLayer);
        _trajectoryType->addAttValue("DrawAs", G4String("Line"));

        _trajectoryType->addAttValue("LineWidthMultiplier", 2.0);
        
        // attributes to draw the points of a track as markers.
        _trajectoryType->addAttValue("MarkName", G4String("Box"));
        _trajectoryType->addAttValue("MarkSize", 4);
        _trajectoryType->addAttValue("MarkType", G4String("Symbol"));
        _trajectoryType->addAttValue("Fill", true);
    }
    return _trajectoryType;
}

HepRepType* G4HepRepSceneHandler::getHitType() {
    if (_hitType == NULL) {
        _hitType = factory->createHepRepType(getEventType(), "Hit");
        _hitType->addAttValue("Layer", hitLayer);
        _hitType->addAttValue("DrawAs", G4String("Point"));
        _hitType->addAttValue("MarkName", G4String("Box"));
        _hitType->addAttValue("MarkSize", 4.0);
        _hitType->addAttValue("MarkType", G4String("Symbol"));
        _hitType->addAttValue("Fill", true);
    }
    return _hitType;
}

HepRepType* G4HepRepSceneHandler::getCalHitType() {
    if (_calHitType == NULL) {
        _calHitType = factory->createHepRepType(getEventType(), "CalHit");
        _calHitType->addAttValue("Layer", calHitLayer);
        _calHitType->addAttValue("Fill", true);
        _calHitType->addAttValue("DrawAs", G4String("Polygon"));
    }
    return _calHitType;
}

HepRepType* G4HepRepSceneHandler::getCalHitFaceType() {
    if (_calHitFaceType == NULL) {
        _calHitFaceType = factory->createHepRepType(getCalHitType(), "CalHitFace");
        _calHitFaceType->addAttValue("PickParent", true);
    }
    return _calHitFaceType;
}

