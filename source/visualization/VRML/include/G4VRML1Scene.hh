// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML1Scene.hh,v 1.1 1999-01-07 16:15:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML1Scene.hh
// Yasuhide Sawada & Satoshi Tanaka

#ifdef  G4VIS_BUILD_VRML_DRIVER

#ifndef G4VRML1_SCENE_HH
#define G4VRML1_SCENE_HH

#include <rw/tvhdict.h>

#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4VScene.hh"

#include "G4FRClient.hh"

class G4VRML1;
class G4VisAttributes;

class G4VRML1Scene: public G4VScene {

	enum { MAX_CONNECTION_TRIAL = 10 } ;

// methods (public) 
public:
	G4VRML1Scene(G4VRML1& system, const G4String& name = "");
	~G4VRML1Scene();
	void AddThis(const G4Box&);
	void AddThis(const G4Cons&);
	void AddThis(const G4Tubs&);
	void AddThis(const G4Trd&);
	void AddThis(const G4Trap&);
	void AddThis(const G4Sphere&);
	void AddThis(const G4Para&);
	void AddThis(const G4Torus&);
	void AddThis(const G4VSolid&);

	void BeginPrimitives(const G4Transform3D& objectTransformation);
	void EndPrimitives();

	void AddPrimitive(const G4Polyline&);
	void AddPrimitive(const G4Polyhedron&);
	void AddPrimitive(const G4NURBS&); 
	void AddPrimitive(const G4Text&); 
	void AddPrimitive(const G4Circle&);
	void AddPrimitive(const G4Square&);
	void AddPrimitive (const G4Polymarker& polymarker)
		{ G4VScene::AddPrimitive (polymarker); }

	void ClearStore();

	void BeginModeling();
	void EndModeling();

	void beginSending();
	void endSending();

	static G4int GetSceneCount() { return fSceneCount; }

	void connectPort(int max_trial = MAX_CONNECTION_TRIAL );
	void closePort();

// methods (private) 
private:
	void SendMaterialNode( const G4VisAttributes*  pAV ); 
	void SendMaterialNode();
	void SendMatrixTransformNode(const G4Transform3D *);
	void SendCubeNode(G4double, G4double, G4double);
	void SendCylinderNode(G4double, G4double);
	void SendSphereNode(G4double);

	void      SendMarkerColor         ( const G4VMarker&  mark ) ;
	void      SendMarkerWorldPosition ( const G4VMarker&  mark ) ;
	G4double  GetMarkerHalfSize       ( const G4VMarker&  mark ) ;

// data 
private:
	G4String fCurrentDEF;

	G4VRML1& fSystem;	// Graphics system for this scene.
	G4FRClient fDest ;

	static G4int fSceneIdCount;
	static G4int fSceneCount;    // No. of existing scenes.

};

#endif //G4VRML1_SCENE_HH
#endif //G4VIS_BUILD_VRML_DRIVER
