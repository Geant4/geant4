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
// $Id: G4HepRepSceneHandler.cc 104290 2017-05-23 13:24:52Z gcosmo $
//

/**
 * @author Mark Donszelmann
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
#include <cassert>

//HepRep
#include "HEPREP/HepRep.h"
#include "G4HepRepMessenger.hh"

//G4
#include "G4PhysicalConstants.hh"
#include "G4Vector3D.hh"
#include "G4Version.hh"
#include "G4Types.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4Polyhedron.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Text.hh"
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
#include "G4AttCheck.hh"

// CHepRep
#include "cheprep/XMLHepRepFactory.h"

// This
#include "G4HepRep.hh"
#include "G4HepRepSceneHandler.hh"
#include "G4HepRepViewer.hh"


using namespace HEPREP;
using namespace cheprep;
using namespace std;

G4int G4HepRepSceneHandler::sceneIdCount = 0;

//#define LDEBUG 1
//#define SDEBUG 1
//#define PDEBUG 1

G4HepRepSceneHandler::G4HepRepSceneHandler (G4VGraphicsSystem& system, const G4String& name)
        : G4VSceneHandler (system, sceneIdCount++, name),
          out                   (0),
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
          writeBinary           (false),
          writeZip              (false),
          writeGZ               (false),
          writeMultipleFiles    (false),
          currentHit            (NULL),
          currentTrack          (NULL),
          _heprep               (NULL),
          _heprepGeometry       (NULL)
{

#ifdef LDEBUG
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
#ifdef LDEBUG
    cout << "G4HepRepSceneHandler::~G4HepRepSceneHandler() " << endl;
#endif
    close();

    delete factory;
    factory = NULL;

    G4HepRep* pHepRepSystem = dynamic_cast<G4HepRep*>(GetGraphicsSystem());
    if (pHepRepSystem) pHepRepSystem->removeSceneHandler();
}


void G4HepRepSceneHandler::open(G4String name) {
    if (writer != NULL) return;

    if (name == "stdout") {
#ifdef LDEBUG
        cout << "G4HepRepSceneHandler::Open() stdout" << endl;
#endif
        writer = factory->createHepRepWriter(&cout, false, false);
        out  = NULL;
        baseName = name;
        eventNumberPrefix = "";
        eventNumberSuffix = "";
        extension = "";
        writeBinary = false;
        writeZip = false;
        writeGZ = false;
        writeMultipleFiles = false;
        eventNumber = 0;
        eventNumberWidth = 0;        
    } else if (name == "stderr") {
#ifdef LDEBUG
        cout << "G4HepRepSceneHandler::Open() stderr" << endl;
#endif
        writer = factory->createHepRepWriter(&cerr, false, false);
        out = NULL;
        baseName = name;
        eventNumberPrefix = "";
        eventNumberSuffix = "";
        extension = "";
        writeBinary = false;
        writeZip = false;
        writeGZ = false;
        writeMultipleFiles = false;
        eventNumber = 0;
        eventNumberWidth = 0;        
    } else {
#ifdef LDEBUG
        cout << "G4HepRepSceneHandler::Open() " << name << endl;
#endif
        if (eventNumberWidth < 0) {
            // derive filename(s)
            // check for extensions
            const unsigned int numberOfExtensions = 8;
            string ext[numberOfExtensions] = {".heprep", ".heprep.xml", ".heprep.zip", ".heprep.gz",
                                              ".bheprep", ".bheprep.xml", ".bheprep.zip", ".bheprep.gz"};
            unsigned int i=0;
            while (i < numberOfExtensions) {
                int dot = name.size() - ext[i].size();
                if ((dot >= 0) && 
                    (name.substr(dot, ext[i].size()) == ext[i])) break;
                i++;
            }
        
            if (i != numberOfExtensions) {
                extension = ext[i];
                writeBinary = i >= (numberOfExtensions/2);
                writeZip = (i == 2) || (i == 6);
                writeGZ = (i == 3) || (i == 7);

                int dot = name.length() - extension.length();
                baseName = (dot >= 0) ? name.substr(0, dot) : "";

            } else {
                // Default for no extension  
                extension = ".heprep.zip";
                writeBinary = false;
                writeZip = true;
                writeGZ = false;
                baseName = name;
            }
        
            writeMultipleFiles = false;
            int startDigit = -1; int endDigit = -1;

			G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();

            string suffix = messenger->getEventNumberSuffix();
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
#ifdef LDEBUG
    cout << "G4HepRepSceneHandler::OpenHepRep() " << endl;
#endif

    if (_heprep != NULL) return;

    // all done on demand, once pointers are set to NULL
    _heprepGeometry       = NULL;
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

#ifdef LDEBUG
    cout << "G4HepRepSceneHandler::CloseHepRep() start" << endl;
#endif

    // if this is the final close, then there should not be any event pending to be written.
    if (final) {
        if (_eventInstanceTree != NULL) {
            cerr << "WARNING: you probably used '/vis/viewer/endOfEventAction accumulate' and "
                 << "forgot to call /vis/viewer/update before exit. No event written." << endl;
        }
    } else {
		
		G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();

        // add geometry to the heprep if there is an event (separate geometries are written
        // using DrawView() called from /vis/viewer/flush)
        if (_eventInstanceTree != NULL) {
            GetCurrentViewer()->DrawView();
        
            // couple geometry

            if ( messenger->appendGeometry()) {
                // couple geometry to event if geometry was written
                if ((_geometryInstanceTree != NULL)) {
                    getEventInstanceTree()->addInstanceTree(getGeometryInstanceTree());
                }
            } else {
                char name[128];
                if (writeMultipleFiles) {
                    sprintf(name, "%s%s%s#%s", baseName.c_str(), "-geometry", extension.c_str(), "G4GeometryData");
                } else {
                    sprintf(name, "%s%s#%s", "geometry", (writeBinary ? ".bheprep" : ".heprep"), "G4GeometryData");
                }   
                getEventInstanceTree()->addInstanceTree(factory->createHepRepTreeID(name, "1.0"));
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
        writeLayers(_heprepGeometry);
        writeLayers(_heprep);

        // open heprep file
        if (writer == NULL) {
            open((GetScene() == NULL) ? G4String("G4HepRepOutput.heprep.zip") : GetScene()->GetName());
        }

        // write out separate geometry
        if (! messenger->appendGeometry() && (_heprepGeometry != NULL)) {
            if (writeMultipleFiles) {
                char fileName[128];
                sprintf(fileName, "%s%s%s", baseName.c_str(), "-geometry", extension.c_str());
                openFile(G4String(fileName));
            }

            char name[128];
            sprintf(name, "%s%s", "geometry", (writeBinary ? ".bheprep" : ".heprep"));
            if (!writeMultipleFiles) {
                writer->addProperty("RecordLoop.ignore", name);
            }

            writer->write(_heprepGeometry, G4String(name));
            
            delete _heprepGeometry;
            _heprepGeometry = NULL;
            
            if (writeMultipleFiles) closeFile();
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
//        eventName << "event-" << setw(eventNumberWidth) << setfill('0') << eventNumber << (writeBinary ? ".bheprep" : ".heprep");
//        writer->write(_heprep, eventName.str());
// Use instead:
        char eventName[128];
        char eventFormat[128];
        sprintf(eventFormat, "%s%d%s%s", "event-%0", eventNumberWidth, "d", (writeBinary ? ".bheprep" : ".heprep"));
        sprintf(eventName, eventFormat, eventNumber);
        if (writer) writer->write(_heprep, G4String(eventName));

        eventNumber++;
    }
    
    delete _heprep;
    _heprep = NULL;

    if (writeMultipleFiles) closeFile();

    return true;
}


void G4HepRepSceneHandler::close() {

#ifdef LDEBUG
    cout << "G4HepRepSceneHandler::Close() " << endl;
#endif

    if (writer == NULL) return;
    
    if (!writeMultipleFiles) {
        closeHepRep(true);
        closeFile();
    }
    
    G4HepRepViewer* viewer = dynamic_cast<G4HepRepViewer*>(GetCurrentViewer());
    viewer->reset();
}

void G4HepRepSceneHandler::openFile(G4String name) {
    out = new ofstream(name.c_str(), std::ios::out | std::ios::binary );
    writer = factory->createHepRepWriter(out, writeZip, writeZip || writeGZ);
}

void G4HepRepSceneHandler::closeFile() {
    writer->close();
    delete writer;
    writer = NULL;

    delete out;
    out = NULL;
}

void G4HepRepSceneHandler::writeLayers(HepRep* heprep) {
    if (heprep == NULL) return;
    heprep->addLayer(geometryLayer); 
    heprep->addLayer(eventLayer);
    heprep->addLayer(calHitLayer);
    heprep->addLayer(trajectoryLayer);
    heprep->addLayer(hitLayer);
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

void G4HepRepSceneHandler::AddSolid(const G4Box& box) {
#ifdef SDEBUG
    cout << "G4HepRepSceneHandler::AddSolid(const G4Box& box)" << endl;
#endif

    if (dontWrite()) return;
	
	G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();

    if (! messenger->useSolids()) {
        G4VSceneHandler::AddSolid(box);
        return;
    }
    
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
    
    vertex1 = (transform) * vertex1;
    vertex2 = (transform) * vertex2;
    vertex3 = (transform) * vertex3;
    vertex4 = (transform) * vertex4;
    vertex5 = (transform) * vertex5;
    vertex6 = (transform) * vertex6;
    vertex7 = (transform) * vertex7;
    vertex8 = (transform) * vertex8;

    HepRepInstance* instance = getGeometryOrEventInstance(getCalHitType());
    addAttributes(instance, getCalHitType());

    setAttribute(instance, "DrawAs", G4String("Prism"));
        
    setVisibility(instance, box);
    setLine(instance, box);
    setColor(instance, getColorFor(box));

    factory->createHepRepPoint(instance, vertex1.x(), vertex1.y(), vertex1.z());
    factory->createHepRepPoint(instance, vertex2.x(), vertex2.y(), vertex2.z());
    factory->createHepRepPoint(instance, vertex3.x(), vertex3.y(), vertex3.z());
    factory->createHepRepPoint(instance, vertex4.x(), vertex4.y(), vertex4.z());
    factory->createHepRepPoint(instance, vertex5.x(), vertex5.y(), vertex5.z());
    factory->createHepRepPoint(instance, vertex6.x(), vertex6.y(), vertex6.z());
    factory->createHepRepPoint(instance, vertex7.x(), vertex7.y(), vertex7.z());
    factory->createHepRepPoint(instance, vertex8.x(), vertex8.y(), vertex8.z());
}


void G4HepRepSceneHandler::AddSolid(const G4Cons& cons) {
#ifdef SDEBUG
    cout << "G4HepRepSceneHandler::AddSolid(const G4Cons& cons)" << endl;
#endif

    if (dontWrite()) return;
	
	G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();

    if (! messenger->useSolids() || (cons.GetDeltaPhiAngle() < twopi)) {
        G4VSceneHandler::AddSolid(cons);
        return;
    }

    G4PhysicalVolumeModel* pPVModel =
      dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
    if (!pPVModel) {
      G4VSceneHandler::AddSolid(cons);
        return;
    }

    G4LogicalVolume* pCurrentLV = pPVModel->GetCurrentLV();
    G4int currentDepth = pPVModel->GetCurrentDepth();
    G4Material* pCurrentMaterial = pPVModel->GetCurrentMaterial();

    G4Point3D vertex1(G4Point3D( 0., 0., cons.GetZHalfLength()));
    G4Point3D vertex2(G4Point3D( 0., 0.,-cons.GetZHalfLength()));

    vertex1 = (transform) * vertex1;
    vertex2 = (transform) * vertex2;

    HepRepInstance* instance = getGeometryInstance(pCurrentLV, pCurrentMaterial, currentDepth);
    setAttribute(instance, "DrawAs", G4String("Cylinder"));
        
    setVisibility(instance, cons);
    setLine(instance, cons);
    setColor(instance, getColorFor(cons));

    HepRepType* type = getGeometryType(pCurrentLV->GetName(), currentDepth);

    // Outer cylinder.
    HepRepInstance* outer = factory->createHepRepInstance(instance, type);
    outer->addAttValue("pickParent",true);
    outer->addAttValue("showParentAttributes",true);
    
    HepRepPoint* op1 = factory->createHepRepPoint(outer, vertex1.x(), vertex1.y(), vertex1.z());
    op1->addAttValue("Radius",cons.GetOuterRadiusPlusZ());
    
    HepRepPoint* op2 = factory->createHepRepPoint(outer, vertex2.x(), vertex2.y(), vertex2.z());
    op2->addAttValue("Radius",cons.GetOuterRadiusMinusZ());

    // Inner cylinder.
    HepRepInstance* inner = factory->createHepRepInstance(instance, type);
    inner->addAttValue("pickParent",true);
    inner->addAttValue("showParentAttributes",true);
    
    HepRepPoint* ip1 = factory->createHepRepPoint(inner, vertex1.x(), vertex1.y(), vertex1.z());
    ip1->addAttValue("Radius",cons.GetInnerRadiusPlusZ());
    
    HepRepPoint* ip2 = factory->createHepRepPoint(inner, vertex2.x(), vertex2.y(), vertex2.z());
    ip2->addAttValue("Radius",cons.GetInnerRadiusMinusZ());
}


void G4HepRepSceneHandler::AddSolid(const G4Tubs& tubs) {
#ifdef SDEBUG
    cout << "G4HepRepSceneHandler::AddSolid(const G4Tubs& tubs)" << endl;
#endif

    if (dontWrite()) return;
	
	G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();

    if (! messenger->useSolids() || (tubs.GetDeltaPhiAngle() < twopi)) {
        G4VSceneHandler::AddSolid(tubs);
        return;
    }
    
    G4PhysicalVolumeModel* pPVModel =
      dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
    if (!pPVModel) {
      G4VSceneHandler::AddSolid(tubs);
        return;
    }

    G4LogicalVolume* pCurrentLV = pPVModel->GetCurrentLV();
    G4int currentDepth = pPVModel->GetCurrentDepth();
    G4Material* pCurrentMaterial = pPVModel->GetCurrentMaterial();

    G4Point3D vertex1(G4Point3D( 0., 0., tubs.GetZHalfLength()));
    G4Point3D vertex2(G4Point3D( 0., 0.,-tubs.GetZHalfLength()));

    vertex1 = (transform) * vertex1;
    vertex2 = (transform) * vertex2;

    HepRepInstance* instance = getGeometryInstance(pCurrentLV, pCurrentMaterial, currentDepth);
    setAttribute(instance, "DrawAs", G4String("Cylinder"));
        
    setVisibility(instance, tubs);
    setLine(instance, tubs);
    setColor(instance, getColorFor(tubs));

    HepRepType* type = getGeometryType(pCurrentLV->GetName(), currentDepth);

    // Outer cylinder.
    HepRepInstance* outer = factory->createHepRepInstance(instance, type);
    outer->addAttValue("Radius",tubs.GetOuterRadius());
    outer->addAttValue("pickParent",true);
    outer->addAttValue("showParentAttributes",true);
    factory->createHepRepPoint(outer, vertex1.x(), vertex1.y(), vertex1.z());
    factory->createHepRepPoint(outer, vertex2.x(), vertex2.y(), vertex2.z());

    // Inner cylinder.
    if (tubs.GetInnerRadius() > 0.) {
        HepRepInstance* inner = factory->createHepRepInstance(instance, type);
        inner->addAttValue("Radius",tubs.GetInnerRadius());
        inner->addAttValue("pickParent",true);
        inner->addAttValue("showParentAttributes",true);
        factory->createHepRepPoint(inner, vertex1.x(), vertex1.y(), vertex1.z());
        factory->createHepRepPoint(inner, vertex2.x(), vertex2.y(), vertex2.z());
    }
}


void G4HepRepSceneHandler::AddSolid(const G4Trd& trd) {
#ifdef SDEBUG
    cout << "G4HepRepSceneHandler::AddSolid(const G4Trd& trd)" << endl;
#endif
    if (dontWrite()) return;
	
	G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();

    if (! messenger->useSolids()) {
        G4VSceneHandler::AddSolid(trd);
        return;
    }
    
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
    
    vertex1 = (transform) * vertex1;
    vertex2 = (transform) * vertex2;
    vertex3 = (transform) * vertex3;
    vertex4 = (transform) * vertex4;
    vertex5 = (transform) * vertex5;
    vertex6 = (transform) * vertex6;
    vertex7 = (transform) * vertex7;
    vertex8 = (transform) * vertex8;

    HepRepInstance* instance = getGeometryOrEventInstance(getCalHitType());

    addAttributes(instance, getCalHitType());

    setAttribute(instance, "DrawAs", G4String("Prism"));
        
    setVisibility(instance, trd);
    setLine(instance, trd);
    setColor(instance, getColorFor(trd));

    factory->createHepRepPoint(instance, vertex1.x(), vertex1.y(), vertex1.z());
    factory->createHepRepPoint(instance, vertex2.x(), vertex2.y(), vertex2.z());
    factory->createHepRepPoint(instance, vertex3.x(), vertex3.y(), vertex3.z());
    factory->createHepRepPoint(instance, vertex4.x(), vertex4.y(), vertex4.z());
    factory->createHepRepPoint(instance, vertex5.x(), vertex5.y(), vertex5.z());
    factory->createHepRepPoint(instance, vertex6.x(), vertex6.y(), vertex6.z());
    factory->createHepRepPoint(instance, vertex7.x(), vertex7.y(), vertex7.z());
    factory->createHepRepPoint(instance, vertex8.x(), vertex8.y(), vertex8.z());
}

void G4HepRepSceneHandler::AddSolid (const G4Trap& trap) { 
    if (dontWrite()) return;
    G4VSceneHandler::AddSolid (trap); 
}

void G4HepRepSceneHandler::AddSolid (const G4Sphere& sphere) { 
    if (dontWrite()) return;
    G4VSceneHandler::AddSolid (sphere); 
}

void G4HepRepSceneHandler::AddSolid (const G4Para& para) { 
    if (dontWrite()) return;
    G4VSceneHandler::AddSolid (para); 
}

void G4HepRepSceneHandler::AddSolid (const G4Torus& torus) { 
    if (dontWrite()) return;
    G4VSceneHandler::AddSolid (torus); 
}

void G4HepRepSceneHandler::AddSolid (const G4Polycone& polycone) {
  if (dontWrite()) return;
  G4VSceneHandler::AddSolid (polycone);
}

void G4HepRepSceneHandler::AddSolid (const G4Polyhedra& polyhedra) {
  if (dontWrite()) return;
  G4VSceneHandler::AddSolid (polyhedra);
}

void G4HepRepSceneHandler::AddSolid (const G4Orb& orb) {
  if (dontWrite()) return;
  G4VSceneHandler::AddSolid (orb);
}

void G4HepRepSceneHandler::AddSolid (const G4Ellipsoid& ellipsoid) {
  if (dontWrite()) return;
  G4VSceneHandler::AddSolid (ellipsoid);
}

void G4HepRepSceneHandler::AddSolid (const G4VSolid& solid) {
    if (dontWrite()) return;
    G4VSceneHandler::AddSolid(solid); 
}


void G4HepRepSceneHandler::AddPrimitive (const G4Polyline& line) {

#ifdef PDEBUG
    cout << "G4HepRepSceneHandler::AddPrimitive(G4Polyline&) " << line.size() << endl;
#endif
    if (dontWrite()) return;

    if (fProcessing2D) {
      static G4bool warned = false;
      if (!warned) {
	warned = true;
	G4Exception
	  ("G4HepRepSceneHandler::AddPrimitive (const G4Polyline&)",
	   "vis-HepRep1001", JustWarning,
	   "2D polylines not implemented.  Ignored.");
      }
      return;
    }

    HepRepInstance* instance = factory->createHepRepInstance(getEventInstance(), getTrajectoryType());

    addAttributes(instance, getTrajectoryType());

    fpVisAttribs = line.GetVisAttributes();
    setColor(instance, GetColor());

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
    if (dontWrite()) return;

    if (fProcessing2D) {
      static G4bool warned = false;
      if (!warned) {
	warned = true;
	G4Exception
	  ("G4HepRepSceneHandler::AddPrimitive (const G4Polymarker&)",
	   "vis-HepRep1002", JustWarning,
	   "2D polymarkers not implemented.  Ignored.");
      }
      return;
    }

    HepRepInstance* instance = factory->createHepRepInstance(getEventInstance(), getHitType());

    addAttributes(instance, getHitType());

    fpVisAttribs = line.GetVisAttributes();
    setColor(instance, GetColor());

    setVisibility(instance, line);

    setMarker(instance, line);

    // Default MarkName is set to Circle for this Type.
    int mtype = line.GetMarkerType();
    
    // Cannot be case statement since line.xxx is not a constant
    if (mtype == line.dots) {
        setAttribute(instance, "Fill", true);
        setColor(instance, GetColor(), G4String("FillColor"));
    } else if (mtype == line.circles) {
    } else if (mtype == line.squares) {
        setAttribute(instance, "MarkName", G4String("Box"));
    } else {
        // line.line + default
        setAttribute(instance, "MarkName", G4String("Plus"));
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
    if (dontWrite()) return;

    if (fProcessing2D) {
      static G4bool warned = false;
      if (!warned) {
	warned = true;
	G4Exception
	  ("G4HepRepSceneHandler::AddPrimitive (const G4Circle&)",
	   "vis-HepRep1003", JustWarning,
	   "2D circles not implemented.  Ignored.");
      }
      return;
    }

    HepRepInstance* instance = factory->createHepRepInstance(getEventInstance(), getHitType());

    addAttributes(instance, getHitType());

    G4Point3D center = transform * circle.GetPosition();

    fpVisAttribs = circle.GetVisAttributes();
    setColor (instance, GetColor());

    setVisibility(instance, circle);

    setMarker(instance, circle);

    factory->createHepRepPoint(instance, center.x(), center.y(), center.z());
}


void G4HepRepSceneHandler::AddPrimitive (const G4Polyhedron& polyhedron) {

#ifdef PDEBUG
    cout << "G4HepRepSceneHandler::AddPrimitive(G4Polyhedron&) " << endl;
#endif
    if (dontWrite()) return;

    if (fProcessing2D) {
      static G4bool warned = false;
      if (!warned) {
	warned = true;
	G4Exception
	  ("G4HepRepSceneHandler::AddPrimitive (const G4Polyhedron&)",
	   "vis-HepRep1004", JustWarning,
	   "2D polyhedra not implemented.  Ignored.");
      }
      return;
    }

    G4Normal3D surfaceNormal;
    G4Point3D vertex;

    if (polyhedron.GetNoFacets()==0) return;

    HepRepInstance* instance = getGeometryOrEventInstance(getCalHitType());
        
    addAttributes(instance, getCalHitType());

    setVisibility(instance, polyhedron);
	
    G4int currentDepth = 0;
    G4PhysicalVolumeModel* pPVModel =
      dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
    if (pPVModel) currentDepth = pPVModel->GetCurrentDepth();

    G4bool notLastFace;
    do {
        HepRepInstance* face;
        if (isEventData()) {
            face = factory->createHepRepInstance(instance, getCalHitFaceType());
        } else {
            face = getGeometryInstance("*Face", currentDepth+1);
            setAttribute(face, "PickParent", true);
            setAttribute(face, "DrawAs", G4String("Polygon"));
        }
        
        setLine(face, polyhedron);
        fpVisAttribs = polyhedron.GetVisAttributes();
        setColor(face, GetColor());
        if (isEventData()) setColor(face, GetColor(), G4String("FillColor"));

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
    if (dontWrite()) return;

    /*** You may need this
    if (fProcessing2D) {
      static G4bool warned = false;
      if (!warned) {
	warned = true;
	G4Exception
	  ("G4HepRepSceneHandler::AddPrimitive (const G4Text&)",
	   "vis-HepRep1005", JustWarning,
	   "2D text not implemented.  Ignored.");
      }
      return;
    }
    ***/

    cout << "G4HepRepSceneHandler::AddPrimitive G4Text : not yet implemented. " << endl;
}


void G4HepRepSceneHandler::AddPrimitive (const G4Square& square) {
#ifdef PDEBUG
    cout << "G4HepRepSceneHandler::AddPrimitive(G4Square&) " << endl;
#endif
    if (dontWrite()) return;

    if (fProcessing2D) {
      static G4bool warned = false;
      if (!warned) {
	warned = true;
	G4Exception
	  ("G4HepRepSceneHandler::AddPrimitive (const G4Square&)",
	   "vis-HepRep1006", JustWarning,
	   "2D squares not implemented.  Ignored.");
      }
      return;
    }

    HepRepInstance* instance = factory->createHepRepInstance(getEventInstance(), getHitType());

    addAttributes(instance, getHitType());

    G4Point3D center = transform * square.GetPosition();

    setColor (instance, getColorFor(square));

    setVisibility(instance, square);

    setMarker(instance, square);

    factory->createHepRepPoint(instance, center.x(), center.y(), center.z());
}

void G4HepRepSceneHandler::AddPrimitive (const G4Scale& scale) { 
    if (dontWrite()) return;
    G4VSceneHandler::AddPrimitive(scale); 
}

void G4HepRepSceneHandler::AddCompound (const G4VTrajectory& trajectory) {
#ifdef PDEBUG
    cout << "G4HepRepSceneHandler::AddCompound(G4VTrajectory&) " << endl;
#endif
    if (dontWrite()) return;
    
    currentTrack = &trajectory;
    G4VSceneHandler::AddCompound(trajectory); 
    currentTrack = NULL;    
}


void G4HepRepSceneHandler::AddCompound (const G4VHit& hit) { 
#ifdef PDEBUG
    cout << "G4HepRepSceneHandler::AddCompound(G4VHit&) " << endl;
#endif
    if (dontWrite()) return;
    
    currentHit = &hit;
    G4VSceneHandler::AddCompound(hit); 
    currentHit = NULL;
}

void G4HepRepSceneHandler::PreAddSolid (const G4Transform3D& objectTransformation,
				  const G4VisAttributes& visAttribs) {

    G4VSceneHandler::PreAddSolid (objectTransformation, visAttribs);

    transform = objectTransformation;
#ifdef SDEBUG
    cout << "G4HepRepSceneHandler::PreAddSolid(G4Transform3D&, G4VisAttributes&)" << endl;
#endif
}


void G4HepRepSceneHandler::PostAddSolid () {
#ifdef SDEBUG
    cout << "G4HepRepSceneHandler::PostAddSolid()" << endl;
#endif
    G4VSceneHandler::PostAddSolid();
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


G4bool G4HepRepSceneHandler::dontWrite() {	
	G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();
    return !( messenger->writeInvisibles() || (fpVisAttribs ? (bool)fpVisAttribs->IsVisible() : true));
}

void G4HepRepSceneHandler::setColor (HepRepAttribute *attribute,
				      const G4Color& color,
				      const G4String& key) {
#ifdef CDEBUG
    cout << "G4HepRepSceneHandler::setColor : red : " << color.GetRed ()   <<
                                  " green : " << color.GetGreen () <<
                                  " blue : " << color.GetBlue ()   << endl;
#endif

    setAttribute(attribute, key, color.GetRed(), color.GetGreen(), color.GetBlue(), color.GetAlpha());
}

G4Color G4HepRepSceneHandler::getColorFor (const G4VSolid& /* solid */) {
    // fpVisAttribs has been set for solids.
    return GetColor();
}

G4Color G4HepRepSceneHandler::getColorFor (const G4Visible& visible) {
  fpVisAttribs = visible.GetVisAttributes();
  return GetColor();
}

void G4HepRepSceneHandler::setVisibility (HepRepAttribute *attribute, const G4VSolid& /* solid */) {
    setAttribute(attribute, "Visibility", (fpVisAttribs ? (bool)fpVisAttribs->IsVisible() : true));
}

void G4HepRepSceneHandler::setVisibility ( HepRepAttribute *attribute, const G4Visible& visible) {
    const G4VisAttributes* atts = visible.GetVisAttributes();

    setAttribute(attribute, "Visibility", (atts && (atts->IsVisible()==0)) ? false : true);
}

void G4HepRepSceneHandler::setLine (HepRepAttribute *attribute, const G4VSolid& /* solid*/) {
    setAttribute(attribute, "LineWidth", 1.0);
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
        fpVisAttribs = marker.GetVisAttributes();
        setColor(attribute, GetColor(), G4String("FillColor"));
    }
}

void G4HepRepSceneHandler::addAttributes(HepRepInstance* instance, HepRepType* type) {
    if (currentHit) {
        vector<G4AttValue>* hitAttValues = currentHit->CreateAttValues();
        const map<G4String,G4AttDef>* hitAttDefs = currentHit->GetAttDefs();

        addAttDefs(getHitType(), hitAttDefs);

        // these attValues are non-standard, so can only be added when we have the attDef.
        type->addAttValue("LVol", G4String(""));
        type->addAttValue("HitType", G4String(""));
        type->addAttValue("ID", -1);
        type->addAttValue("Column", -1);
        type->addAttValue("Row", -1);
        type->addAttValue("Energy", 0.0);
        type->addAttValue("Pos", G4String(""));

        addAttVals(instance, hitAttDefs, hitAttValues);
    
        delete hitAttValues;
    
    } else if (currentTrack) {
        vector<G4AttValue>* trajectoryAttValues = currentTrack->CreateAttValues();
        const map<G4String,G4AttDef>* trajectoryAttDefs = currentTrack->GetAttDefs();

        addAttDefs(type, trajectoryAttDefs);
    
        // these attValues are non-standard, so can only be added when we have the attDef.
        type->addAttValue("Ch", 0.0);
        type->addAttValue("Color", 1.0, 1.0, 1.0, 1.0);
        type->addAttValue("ID", -1);
        type->addAttValue("IMom", G4String(""));
        type->addAttValue("IMag", 0.0);
        type->addAttValue("PDG", -1);
        type->addAttValue("PN", G4String(""));
        type->addAttValue("PID", -1);

        addAttVals(instance, trajectoryAttDefs, trajectoryAttValues);
    
        delete trajectoryAttValues;
        
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
        definition->addAttDef(attDefIterator->first, attDefIterator->second.GetDesc(),
                        attDefIterator->second.GetCategory(), attDefIterator->second.GetExtra());
        attDefIterator++;
    }
}

void G4HepRepSceneHandler::addAttVals(HepRepAttribute* attribute, const map<G4String,G4AttDef>* attDefs, vector<G4AttValue>* attValues) {
    if (attValues == NULL) return;

    // Copy the instance's G4AttValues to HepRepAttValues.
    for (vector<G4AttValue>::iterator attValIterator = attValues->begin(); attValIterator != attValues->end(); attValIterator++) {        
        G4String name = attValIterator->GetName();

        HepRepPoint* point = dynamic_cast<HepRepPoint*>(attribute);
        if ((name == "Pos") && (point != NULL)) {
            G4String pos = attValIterator->GetValue();
//            cout << "Pos* " << pos << endl;
            int is = 0;
            int in = 0;
            int im = 0;
            G4String unit;
            for (unsigned int i=0; i<pos.length(); i++) {
                if (pos[i] == ' ') {
                    if (in == 0) {
                        // first coordinate
                        double factor = atof(pos.substr(is, i-is).c_str())/point->getX();
                        im = (int)(std::log10(factor)+((factor < 1) ? -0.5 : 0.5));
//                        cout << factor << ", " << im << endl;
                    } else if (in == 3) {
                        // unit
                        unit = pos.substr(is, i-is);
                        if (unit == G4String("mum")) {
                            im += -6;
                        } else if (unit == G4String("mm")) {
                            im += -3;
                        } else if (unit == G4String("cm")) {
                            im += -2;
                        } else if (unit == G4String("m")) {
                            im += 0;
                        } else if (unit == G4String("km")) {
                            im += 3;
                        } else {
                            cerr << "HepRepSceneHandler: Unrecognized Unit: '" << unit << "'" << endl;
                        }
                    }
                    is = i+1;
                    in++;
                }
            }
            switch(im) {
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
                    cerr << "HepRepSceneHandler: No valid unit found for im: " << im << endl;
                    unit = G4String("*im");                        
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
    G4PhysicalVolumeModel* pPVModel =
      dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
    return !pPVModel || fReadyForTransients || currentHit || currentTrack;
}

void G4HepRepSceneHandler::addTopLevelAttributes(HepRepType* type) {
    
    // Some non-standard attributes
    type->addAttDef(  "Generator", "Generator of the file", "General", "");
    type->addAttValue("Generator", G4String("Geant4"));

    type->addAttDef(  "GeneratorVersion", "Version of the Generator", "General", "");
    G4String versionString = G4Version;
    versionString = versionString.substr(1,versionString.size()-2);
    versionString = " Geant4 version " + versionString + "   " + G4Date;
    type->addAttValue("GeneratorVersion", versionString);
    
    const G4ViewParameters parameters = GetCurrentViewer()->GetViewParameters();
    const G4Vector3D& viewPointDirection = parameters.GetViewpointDirection();    
    type->addAttDef(  "ViewTheta", "Theta of initial suggested viewpoint", "Draw", "rad");
    type->addAttValue("ViewTheta", viewPointDirection.theta());    
    
    type->addAttDef(  "ViewPhi", "Phi of initial suggested viewpoint", "Draw", "rad");
    type->addAttValue("ViewPhi", viewPointDirection.phi());    
    
    type->addAttDef(  "ViewScale", "Scale of initial suggested viewpoint", "Draw", "");
    type->addAttValue("ViewScale", parameters.GetZoomFactor());    
    
// FIXME, no way to set these
    type->addAttDef(  "ViewTranslateX", "Translate in X of initial suggested viewpoint", "Draw", "");
    type->addAttValue("ViewTranslateX", 0.0);    
    
    type->addAttDef(  "ViewTranslateY", "Translate in Y of initial suggested viewpoint", "Draw", "");
    type->addAttValue("ViewTranslateY", 0.0);    
    
    type->addAttDef(  "ViewTranslateZ", "Translate in Z of initial suggested viewpoint", "Draw", "");
    type->addAttValue("ViewTranslateZ", 0.0);    
    
    type->addAttDef(  "PointUnit", "Length", "Physics", "");
    type->addAttValue("PointUnit", G4String("m"));
	
	G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();

    type->addAttDef(  "UseSolids", "Use HepRep Solids rather than Geant4 Primitives", "Draw", "");
    type->addAttValue("UseSolids", messenger->useSolids());

    type->addAttDef(  "WriteInvisibles", "Write Invisible Objects", "Draw", "");
    type->addAttValue("WriteInvisibles", messenger->writeInvisibles());
}           


HepRepInstance* G4HepRepSceneHandler::getGeometryOrEventInstance(HepRepType* type) {
    if (isEventData()) {
      return factory->createHepRepInstance(getEventInstance(), type);
    } else {
      G4PhysicalVolumeModel* pPVModel =
	dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
      assert(pPVModel);  // To keep Coverity happy.
      G4LogicalVolume* pCurrentLV = pPVModel->GetCurrentLV();
      G4int currentDepth = pPVModel->GetCurrentDepth();
      G4Material* pCurrentMaterial = pPVModel->GetCurrentMaterial();
      return getGeometryInstance(pCurrentLV, pCurrentMaterial, currentDepth);
    }
}

HepRep* G4HepRepSceneHandler::getHepRep() {
    if (_heprep == NULL) {
        // Create the HepRep that holds the Trees.
        _heprep = factory->createHepRep();
    }
    return _heprep;
}   

HepRep* G4HepRepSceneHandler::getHepRepGeometry() {
    if (_heprepGeometry == NULL) {
        // Create the HepRep that holds the Trees.
        _heprepGeometry = factory->createHepRep();
    }
    return _heprepGeometry;
}   

HepRepInstanceTree* G4HepRepSceneHandler::getGeometryInstanceTree() {
    if (_geometryInstanceTree == NULL) {
        // Create the Geometry InstanceTree.
        _geometryInstanceTree = factory->createHepRepInstanceTree("G4GeometryData", "1.0", getGeometryTypeTree());
		
		G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();
        if ( messenger->appendGeometry()) {
            getHepRep()->addInstanceTree(_geometryInstanceTree);
        } else {
            getHepRepGeometry()->addInstanceTree(_geometryInstanceTree);
        }
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

HepRepInstance* G4HepRepSceneHandler::getGeometryInstance(G4LogicalVolume* volume, G4Material* material, int depth) {
    HepRepInstance* instance = getGeometryInstance(volume->GetName(), depth);

    setAttribute(instance, "LVol",       volume->GetName());
    G4Region* region = volume->GetRegion();
    G4String regionName = region? region->GetName(): G4String("No region");
    setAttribute(instance, "Region",     regionName);
    setAttribute(instance, "RootRegion", volume->IsRootRegion());
    setAttribute(instance, "Solid",      volume->GetSolid()->GetName());
    setAttribute(instance, "EType",      volume->GetSolid()->GetEntityType());
    G4String matName = material? material->GetName(): G4String("No material");
    setAttribute(instance, "Material",   matName );
    G4double matDensity = material? material->GetDensity(): 0.;
    setAttribute(instance, "Density",    matDensity);
    G4double matRadlen = material? material->GetRadlen(): 0.;
    setAttribute(instance, "Radlen",     matRadlen);
    
    G4State matState = material? material->GetState(): kStateUndefined;
    G4String state = materialState[matState];
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
		
		G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();
        if ( messenger->appendGeometry()) {
            getHepRep()->addTypeTree(_geometryTypeTree);
        } else {
            getHepRepGeometry()->addTypeTree(_geometryTypeTree);
        }
    }
    return _geometryTypeTree;
}

HepRepType* G4HepRepSceneHandler::getGeometryRootType() {
    if (_geometryRootType == NULL) {
        // Create the top level Geometry Type.
        _geometryRootType = factory->createHepRepType(getGeometryTypeTree(), rootVolumeName);
        _geometryRootType->addAttValue("Layer", geometryLayer);
    
        // Add attdefs used by all geometry types.
        _geometryRootType->addAttDef  ("LVol", "Logical Volume", "Physics","");
        _geometryRootType->addAttValue("LVol", G4String(""));
        _geometryRootType->addAttDef  ("Region", "Cuts Region", "Physics","");
        _geometryRootType->addAttValue("Region", G4String(""));
        _geometryRootType->addAttDef  ("RootRegion", "Root Region", "Physics","");
        _geometryRootType->addAttValue("RootRegion", false);
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

        addTopLevelAttributes(_geometryRootType);       
                
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
        _eventType = factory->createHepRepType(getEventTypeTree(), "Event");
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

        addTopLevelAttributes(_eventType);
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

