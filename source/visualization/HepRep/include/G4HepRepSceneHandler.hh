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
 */

#ifndef G4HEPREPSCENEHANDLER_HH
#define G4HEPREPSCENEHANDLER_HH 1

#include "globals.hh"
#include <iostream>
#include <stack>
#include <map>

// HepRep
#include "HEPREP/HepRep.h"

//G4
#include "G4VGraphicsSystem.hh"
#include "G4VSceneHandler.hh"
#include "G4Visible.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalVolumeModel.hh"


//class G4HepRep;

class G4HepRepSceneHandler: public G4VSceneHandler {

    public:
        G4HepRepSceneHandler (G4VGraphicsSystem& system, const G4String& name = "");
        virtual ~G4HepRepSceneHandler ();

        void AddThis (const G4Box& box)                 { G4VSceneHandler::AddThis (box); }
        void AddThis (const G4Cons& cons)               { G4VSceneHandler::AddThis (cons); }
        void AddThis (const G4Tubs& tubs)               { G4VSceneHandler::AddThis (tubs); }
        void AddThis (const G4Trd& trd)                 { G4VSceneHandler::AddThis (trd); }
        void AddThis (const G4Trap& trap)               { G4VSceneHandler::AddThis (trap); }
        void AddThis (const G4Sphere& sphere)           { G4VSceneHandler::AddThis (sphere); }
        void AddThis (const G4Para& para)               { G4VSceneHandler::AddThis (para); }
        void AddThis (const G4Torus& torus)             { G4VSceneHandler::AddThis (torus); }
        void AddThis (const G4Polycone& polycone)       { G4VSceneHandler::AddThis (polycone); }
        void AddThis (const G4Polyhedra& polyhedra)     { G4VSceneHandler::AddThis (polyhedra); }
        void AddThis (const G4VSolid& solid)            { G4VSceneHandler::AddThis(solid); }
        void AddThis (const G4VTrajectory&);
        void AddThis (const G4VHit& hit)                { G4VSceneHandler::AddThis(hit); }

        void PreAddThis (const G4Transform3D& objectTransformation, const G4VisAttributes& visAttribs);
        void PostAddThis ();

        void EstablishSpecials(G4PhysicalVolumeModel&);

        void AddPrimitive (const G4Polyline&);
        void AddPrimitive (const G4Text&);
        void AddPrimitive (const G4Circle&);
        void AddPrimitive (const G4Square&);
        void AddPrimitive (const G4Polyhedron&);
        void AddPrimitive (const G4NURBS&);

        void AddPrimitive (const G4Polymarker&);
        void AddPrimitive (const G4Scale& scale)            { G4VSceneHandler::AddPrimitive(scale); }

        void BeginPrimitives (const G4Transform3D& objectTransformation);
        void EndPrimitives ();
        void BeginModeling ();
        void EndModeling ();

        static G4int GetSceneCount ()                       { return sceneCount; }

        void OpenHepRep();
        bool CloseHepRep();

    private:
        static G4int sceneCount;
        G4int eventNumber;

        G4int currentDepth;
        G4VPhysicalVolume* currentPV;
        G4LogicalVolume* currentLV;
        const G4ModelingParameters* originalMP;  // Keeps pointer to original.
        G4ModelingParameters* nonCullingMP;      // For temporary non-culling.

        G4Transform3D transform;

//        G4int geomParentDepth;
//        std::map<G4String, HEPREP::HepRepType *> geometryTypeByPath;
//        std::vector<HEPREP::HepRepType *> geometryTypeByDepth;
//        std::stack<HEPREP::HepRepInstance *> geomParentInstanceS;

//        G4int eventParentDepth;
//        std::map<G4String, HEPREP::HepRepType *> eventTypeFullNameMap;
//        std::stack<G4String> eventParentTypeFullNameS;
//        std::stack<HEPREP::HepRepInstance *> eventParentInstanceS;

        std::ostream* out;
        HEPREP::HepRepFactory* factory;
        HEPREP::HepRepWriter* writer;
        HEPREP::HepRep* heprep;

        HEPREP::HepRepType* detectorType;
        HEPREP::HepRepInstanceTree* geometryInstanceTree;
        HEPREP::HepRepType* eventType;
        HEPREP::HepRepInstanceTree* eventInstanceTree;

//        char geomTypeFname [256];
//        char geomInstanceFname [256];
//        char eventTypeFname [256];
//        char eventInstanceFname [256];

        void SetColour (HEPREP::HepRepAttribute *attribute, const G4Colour& color,
			            const G4String& key = G4String("Color"));
        void SetLine (HEPREP::HepRepInstance *instance, const G4Visible& visible);
        void SetMarker (HEPREP::HepRepInstance *instance, const G4VMarker& marker);

        HEPREP::HepRepInstance* CreateGeometryInstance(G4String typeName, G4int depth);
        HEPREP::HepRepInstance* CreateEventInstance(G4String typeName, G4int depth,
					    const std::map<G4String,G4AttDef>* attDefs = NULL,
					    std::vector<G4AttValue>* attValues = NULL);

        bool IsEventData ();

        void Open(G4String name);
        void Close();
};

#endif

