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
#include "g4std/iostream"

// HepRep
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

        HEPREP::HepRep *GetHepRep() { return heprep; }
        HEPREP::HepRepWriter *GetHepRepWriter() { return writer; }
        HEPREP::HepRepFactory *GetHepRepFactory();
        void open();
        void close();

    protected:
//        void RequestPrimitives (const G4VSolid& solid);

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

        void SetColour (HEPREP::HepRepAttribute *attribute, const G4Colour& color,const G4String& key = G4String("Color"));
        void SetLine (HEPREP::HepRepInstance *instance, const G4Visible& visible);
        void SetMarker (HEPREP::HepRepInstance *instance, const G4VMarker& marker);
        HEPREP::HepRepInstance* CreateInstance(HEPREP::HepRepInstance* p, HEPREP::HepRepType* altType);
        bool IsEventData ();

        G4std::ostream* out;
        HEPREP::HepRepFactory* heprepFactory;
        HEPREP::HepRepWriter* writer;
        HEPREP::HepRepInstance *parent;
        HEPREP::HepRep *heprep;
        HEPREP::HepRepType *geometryType;
        HEPREP::HepRepType *eventType;
        HEPREP::HepRepType *trackType;
        HEPREP::HepRepType *calHitType;
        HEPREP::HepRepType *hitType;
};

#endif

