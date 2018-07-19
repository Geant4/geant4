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
// $Id: G4HepRepSceneHandler.hh 99152 2016-09-07 08:04:30Z gcosmo $
//

/**
 * @author Mark Donszelmann
 */

#ifndef G4HEPREPSCENEHANDLER_HH
#define G4HEPREPSCENEHANDLER_HH 1

#include "globals.hh"
#include <iostream>
#include <stack>
#include <map>
#include <vector>

// HepRep
#include "HEPREP/HepRep.h"

//G4
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Sphere.hh"
#include "G4Para.hh"
#include "G4Torus.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"

#include "G4VGraphicsSystem.hh"
#include "G4VSceneHandler.hh"
#include "G4Visible.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalVolumeModel.hh"

class G4HepRepSceneHandler: public G4VSceneHandler {

    public:
        G4HepRepSceneHandler (G4VGraphicsSystem& system, const G4String& name = "");
        virtual ~G4HepRepSceneHandler ();

        void AddSolid (const G4Box& box);
        void AddSolid (const G4Cons& cons);
        void AddSolid (const G4Tubs& tubs);
        void AddSolid (const G4Trd& trd);
        void AddSolid (const G4Trap& trap);          
        void AddSolid (const G4Sphere& sphere);      
        void AddSolid (const G4Para& para);          
        void AddSolid (const G4Torus& torus);        
        void AddSolid (const G4Polycone& polycone);
        void AddSolid (const G4Polyhedra& polyhedra);
        void AddSolid (const G4Orb& orb);
        void AddSolid (const G4Ellipsoid& ellipsoid);
        void AddSolid (const G4VSolid& solid);
        
        void AddCompound (const G4VTrajectory&);
        void AddCompound (const G4VHit& hit);
        void AddCompound (const G4VDigi& digi) {
	  G4VSceneHandler::AddCompound(digi);
	}
        void AddCompound (const G4THitsMap<G4double>& hits) {
	  G4VSceneHandler::AddCompound(hits);
	}
        void AddCompound (const G4THitsMap<G4StatDouble>& hits) {
	  G4VSceneHandler::AddCompound(hits);
	}

        void PreAddSolid (const G4Transform3D& objectTransformation, const G4VisAttributes& visAttribs);
        void PostAddSolid ();

        void AddPrimitive (const G4Polyline&);
        void AddPrimitive (const G4Text&);
        void AddPrimitive (const G4Circle&);
        void AddPrimitive (const G4Square&);
        void AddPrimitive (const G4Polyhedron&);

        void AddPrimitive (const G4Polymarker&);
        void AddPrimitive (const G4Scale& scale);

        void BeginPrimitives (const G4Transform3D& objectTransformation);
        void EndPrimitives ();
        void BeginModeling ();
        void EndModeling ();

        void openHepRep();
        bool closeHepRep(bool final = false);
        void openFile(G4String name);
        void closeFile();

    private:
        G4HepRepSceneHandler (const G4HepRepSceneHandler&);
        G4HepRepSceneHandler& operator= (const G4HepRepSceneHandler&);
        static G4int sceneIdCount;

        G4Transform3D transform;

        std::ostream* out;
        HEPREP::HepRepFactory* factory;
        HEPREP::HepRepWriter* writer;
        
        // Methods
        G4bool dontWrite();
        
        void setColor (HEPREP::HepRepAttribute *attribute, const G4Color& color,
			            const G4String& key = G4String("Color"));
        G4Color getColorFor (const G4Visible& visible);
        G4Color getColorFor (const G4VSolid& solid);
        
        void setVisibility (HEPREP::HepRepAttribute *attribute, const G4VSolid& solid);
        void setLine   (HEPREP::HepRepAttribute *attribute, const G4VSolid& solid);
        
        void setVisibility (HEPREP::HepRepAttribute *attribute, const G4Visible& visible);
        void setLine   (HEPREP::HepRepAttribute *attribute, const G4Visible& visible);

        void setMarker (HEPREP::HepRepAttribute *attribute, const G4VMarker& marker);

        inline void setAttribute (HEPREP::HepRepAttribute* attribute, G4String name, char* value) {
            setAttribute(attribute, name, G4String(value));
        }
        void setAttribute (HEPREP::HepRepAttribute* attribute, G4String name, G4String value);
        void setAttribute (HEPREP::HepRepAttribute* attribute, G4String name, bool value);
        void setAttribute (HEPREP::HepRepAttribute* attribute, G4String name, double value);
        void setAttribute (HEPREP::HepRepAttribute* attribute, G4String name, int value);
        void setAttribute (HEPREP::HepRepAttribute* attribute, G4String name, double red, double green, double blue, double alpha);

        bool isEventData();

        void open(G4String name);
        void close();

        void writeLayers(HEPREP::HepRep* heprep);

        void addAttributes(HEPREP::HepRepInstance* instance, HEPREP::HepRepType* type);

        void addAttDefs(HEPREP::HepRepDefinition* definition, const std::map<G4String,G4AttDef>* attDefs);
        void addAttVals(HEPREP::HepRepAttribute* attribute, const std::map<G4String,G4AttDef>* attDefs, std::vector<G4AttValue>* attValues);

        void addTopLevelAttributes(HEPREP::HepRepType* type);

        HEPREP::HepRepInstance*     getGeometryOrEventInstance(HEPREP::HepRepType* type);

        // Returns the particular instance/type or if not created, creates them and adds them to the HepRep
        HEPREP::HepRep*             getHepRep();
        HEPREP::HepRep*             getHepRepGeometry();
        HEPREP::HepRepInstanceTree* getGeometryInstanceTree();
        HEPREP::HepRepInstance*     getGeometryInstance(G4LogicalVolume* volume, G4Material* material, int depth);
        HEPREP::HepRepInstance*     getGeometryInstance(G4String volumeName, int depth);
        HEPREP::HepRepInstance*     getGeometryRootInstance();
        HEPREP::HepRepTypeTree*     getGeometryTypeTree();
        HEPREP::HepRepType*         getGeometryType(G4String volumeName, int depth);
        HEPREP::HepRepType*         getGeometryRootType();
        HEPREP::HepRepInstanceTree* getEventInstanceTree();
        HEPREP::HepRepInstance*     getEventInstance();
        HEPREP::HepRepTypeTree*     getEventTypeTree();
        HEPREP::HepRepType*         getEventType();
        HEPREP::HepRepType*         getTrajectoryType       ();
        HEPREP::HepRepType*         getHitType              ();
        HEPREP::HepRepType*         getCalHitType           ();
        HEPREP::HepRepType*         getCalHitFaceType       ();

        G4String getFullTypeName(G4String volumeName, int depth);
        G4String getParentTypeName(int currentDepth);

        // initialized Member Variables
        G4String geometryLayer, eventLayer, calHitLayer;
        G4String trajectoryLayer, hitLayer;
        G4String rootVolumeName;
        
        G4String baseName;
        G4String eventNumberPrefix;
        G4String eventNumberSuffix;
        G4int eventNumber;
        G4int eventNumberWidth;
        G4String extension;
        G4bool writeBinary;
        G4bool writeZip;
        G4bool writeGZ;
        G4bool writeMultipleFiles;
        const G4VHit* currentHit;
        const G4VTrajectory* currentTrack;

        // DO NOT USE member vars directly, use get methods.
        HEPREP::HepRep*                         _heprep;
        HEPREP::HepRep*                         _heprepGeometry;
        HEPREP::HepRepInstanceTree*             _geometryInstanceTree;
        std::vector<HEPREP::HepRepInstance*>    _geometryInstance;
        HEPREP::HepRepInstance*                 _geometryRootInstance;
        HEPREP::HepRepTypeTree*                 _geometryTypeTree;
        std::vector<G4String>                   _geometryTypeName;
        std::map<G4String, HEPREP::HepRepType*> _geometryType;
        HEPREP::HepRepType*                     _geometryRootType;
        HEPREP::HepRepInstanceTree*             _eventInstanceTree;
        HEPREP::HepRepInstance*                 _eventInstance;
        HEPREP::HepRepTypeTree*                 _eventTypeTree;
        HEPREP::HepRepType*                     _eventType;
        HEPREP::HepRepType*                     _trajectoryType;
        HEPREP::HepRepType*                     _hitType;
        HEPREP::HepRepType*                     _calHitType;
        HEPREP::HepRepType*                     _calHitFaceType;        

        std::map<int, G4String> materialState;
};

#endif

