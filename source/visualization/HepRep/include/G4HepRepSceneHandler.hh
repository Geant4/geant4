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
#include <vector>

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
        const G4ModelingParameters* originalMP;  // Keeps pointer to original.
        G4ModelingParameters* nonCullingMP;      // For temporary non-culling.

        G4Transform3D transform;

        std::ostream* out;
        HEPREP::HepRepFactory* factory;
        HEPREP::HepRepWriter* writer;
        
        // Methods
        void SetColor (HEPREP::HepRepAttribute *attribute, const G4Color& color,
			            const G4String& key = G4String("Color"));
        void SetLine   (HEPREP::HepRepInstance *instance, const G4Visible& visible);
        void SetMarker (HEPREP::HepRepInstance *instance, const G4VMarker& marker);

        bool IsEventData ();

        void Open(G4String name);
        void Close();

        void addAttDefs(HEPREP::HepRepType* type, const std::map<G4String,G4AttDef>* attDefs);
        void addAttVals(HEPREP::HepRepInstance* instance, std::vector<G4AttValue>* attValues);


        // Returns the particular instance/type or if not created, creates them and adds them to the HepRep
        HEPREP::HepRep*             getHepRep();
        HEPREP::HepRepInstanceTree* getGeometryInstanceTree();
        HEPREP::HepRepInstance*     getGeometryInstance(G4LogicalVolume* volume, int depth);
        HEPREP::HepRepInstance*     getGeometryInstance(G4String volumeName, int depth);
        HEPREP::HepRepInstance*     getGeometryRootInstance();
        HEPREP::HepRepTypeTree*     getGeometryTypeTree();
        HEPREP::HepRepType*         getGeometryType(G4String volumeName, int depth);
        HEPREP::HepRepType*         getGeometryRootType();
        HEPREP::HepRepInstanceTree* getEventInstanceTree();
        HEPREP::HepRepInstance*     getEventInstance();
        HEPREP::HepRepTypeTree*     getEventTypeTree();
        HEPREP::HepRepType*         getEventType();
        HEPREP::HepRepType*         getTrajectoryType       (const std::map<G4String,G4AttDef>* attDefs = NULL);
        HEPREP::HepRepType*         getTrajectoryPointType  (const std::map<G4String,G4AttDef>* attDefs = NULL);
        HEPREP::HepRepType*         getHitType              (const std::map<G4String,G4AttDef>* attDefs = NULL);
        HEPREP::HepRepType*         getCalHitType           (const std::map<G4String,G4AttDef>* attDefs = NULL);
        HEPREP::HepRepType*         getCalHitFaceType       (const std::map<G4String,G4AttDef>* attDefs = NULL);

        G4String getFullTypeName(G4String volumeName, int depth);
        G4String getParentTypeName(int currentDepth);

        // initialized Member Variables
        G4String geometryLayer, eventLayer, calHitLayer;
        G4String trajectoryLayer, trajectoryPointLayer, hitLayer;
        G4String rootVolumeName;
        
        G4int eventNumber;

        G4int currentDepth;
        G4VPhysicalVolume* currentPV;
        G4LogicalVolume* currentLV;
        
        // DO NOT USE member vars directly, use get methods.
        HEPREP::HepRep*                         _heprep;
        HEPREP::HepRepInstanceTree*             _geometryInstanceTree;
        std::stack<HEPREP::HepRepInstance*>     _geometryInstance;
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
        HEPREP::HepRepType*                     _trajectoryPointType;
        HEPREP::HepRepType*                     _hitType;
        HEPREP::HepRepType*                     _calHitType;
        HEPREP::HepRepType*                     _calHitFaceType;        

        std::map<int, G4String> materialState;
};

#endif

