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

#include "XMLHepRepStreamer.h"

//class G4HepRep;

class G4HepRepSceneHandler: public G4VSceneHandler {

    public:
        G4HepRepSceneHandler (G4VGraphicsSystem& system, const G4String& name = "");
        virtual ~G4HepRepSceneHandler ();

        void AddThis (const G4Box&);
        void AddThis (const G4Cons&);
        void AddThis (const G4Tubs&);
        void AddThis (const G4Trd&);
        void AddThis (const G4Trap&);
        void AddThis (const G4Sphere&);
        void AddThis (const G4Para&);
        void AddThis (const G4Torus&);
        void AddThis (const G4Polycone&);
        void AddThis (const G4Polyhedra&);
        void AddThis (const G4VSolid&);
        void AddThis (const G4VTrajectory&);
        void AddThis (const G4VHit&);

        void PreAddThis (const G4Transform3D& objectTransformation,
			             const G4VisAttributes& visAttribs);
        void PostAddThis ();

        void EstablishSpecials(G4PhysicalVolumeModel&);

        void AddPrimitive (const G4Polyline&);
        void AddPrimitive (const G4Text&);
        void AddPrimitive (const G4Circle&);
        void AddPrimitive (const G4Square&);
        void AddPrimitive (const G4Polyhedron&);
        void AddPrimitive (const G4NURBS&);

        void AddPrimitive (const G4Polymarker&);
        void AddPrimitive (const G4Scale& scale) {
            G4VSceneHandler::AddPrimitive(scale);
        }

        void BeginPrimitives (const G4Transform3D& objectTransformation);
        void EndPrimitives ();
        void BeginModeling ();
        void EndModeling ();

        static G4int GetSceneCount () { return sceneCount; }

        void close();
    private:
        static G4int sceneCount;      // No. of extanct scenes.
        static G4int sceneIdCount;    // static counter for scenes.
        G4int fileNo;                   //  file sequence number
        G4int currentDepth;
        G4VPhysicalVolume* currentPV;
        G4LogicalVolume* currentLV;
        const G4ModelingParameters* originalMP;  // Keeps pointer to original.
        G4ModelingParameters* nonCullingMP;      // For temporary non-culling.

        G4Transform3D transform;

        G4int geomParentDepth;
        std::map<G4String, HEPREP::HepRepType *> geomTypeFullNameMap;
        std::stack<G4String> geomParentTypeFullNameS;
        std::stack<HEPREP::HepRepInstance *> geomParentInstanceS;

        G4int eventParentDepth;
        std::map<G4String, HEPREP::HepRepType *> eventTypeFullNameMap;
        std::stack<G4String> eventParentTypeFullNameS;
        std::stack<HEPREP::HepRepInstance *> eventParentInstanceS;

        std::ostream* out;
        HEPREP::HepRepFactory* factory;
        XMLHepRepStreamer* writer;
        HEPREP::HepRep* heprep;

        std::ostream* geomTypeOut;
        std::ostream* geomInstanceOut;
        HEPREP::HepRepFactory* geomTypeHepRepFactory;
        HEPREP::HepRepFactory* geomInstanceHepRepFactory;
        HEPREP::HepRepWriter* geomTypeWriter;
        HEPREP::HepRepWriter* geomInstanceWriter;
        HEPREP::HepRep* geomTypeHeprep;
        HEPREP::HepRep* geomInstanceHeprep;
        HEPREP::HepRepInstanceTree* geomInstanceTree;

        std::ostream* eventTypeOut;
        std::ostream* eventInstanceOut;
        HEPREP::HepRepFactory* eventTypeHepRepFactory;
        HEPREP::HepRepFactory* eventInstanceHepRepFactory;
        HEPREP::HepRepWriter* eventTypeWriter;
        HEPREP::HepRepWriter* eventInstanceWriter;
        HEPREP::HepRep* eventTypeHeprep;
        HEPREP::HepRep* eventInstanceHeprep;
        HEPREP::HepRepInstanceTree* eventInstanceTree;

        char geomTypeFname [256];
        char geomInstanceFname [256];
        char eventTypeFname [256];
        char eventInstanceFname [256];

        void SetColour (HEPREP::HepRepAttribute *attribute, const G4Colour& color,
			            const G4String& key = G4String("Color"));
        void SetLine (HEPREP::HepRepInstance *instance, const G4Visible& visible);
        void SetMarker (HEPREP::HepRepInstance *instance, const G4VMarker& marker);

        HEPREP::HepRepInstance* CreateGeomInstance(G4String typeName, G4int depth);
        HEPREP::HepRepInstance* CreateEventInstance(G4String typeName, G4int depth,
					    const std::map<G4String,G4AttDef>* attDefs = NULL,
					    std::vector<G4AttValue>* attValues = NULL);

        bool IsEventData ();

        void open();
        void openGeomHepRep();
        void openEventHepRep();
        void mergeAndDelete(char* fname);
};

#endif

