// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML2SceneHandler.hh,v 1.4 1999-12-15 14:54:15 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML2SceneHandler.hh
// Satoshi Tanaka & Yasuhide Sawada

#ifdef  G4VIS_BUILD_VRML_DRIVER

#ifndef G4VRML2_SCENE_HANDLER_HH
#define G4VRML2_SCENE_HANDLER_HH

#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4VSceneHandler.hh"

#include "G4FRClient.hh"

class G4VRML2;
class G4VisAttributes;

class G4VRML2SceneHandler: public G4VSceneHandler {

	enum { MAX_CONNECTION_TRIAL = 10 } ;

// methods (public) 
public:
	G4VRML2SceneHandler(G4VRML2& system, const G4String& name = "");
	virtual ~G4VRML2SceneHandler();
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
		{ G4VSceneHandler::AddPrimitive (polymarker); }

	void ClearStore();

	void BeginModeling();
	void EndModeling();

	static G4int GetSceneCount() { return fSceneCount; }

	void VRMLBeginModeling();
	void VRMLEndModeling();

	void connectPort(int max_trial = MAX_CONNECTION_TRIAL );
	void closePort();

// methods (private) 
private:
	void      SendMaterialNode        ( const G4VisAttributes*  pAV ); 
	void      SendMaterialNode        ();

	void      SendLineColor           ( const G4VisAttributes*  pAV ); 
	void      SendMarkerColor         ( const G4VMarker&  mark ) ;
	void      SendMarkerWorldPosition ( const G4VMarker&  mark ) ;

	G4double  GetMarkerHalfSize       ( const G4VMarker&  mark ) ;
	void      GetMarkerWorldPosition  (	const G4VMarker&  mark , 
						double* pX             ,
						double* pY             ,
						double* pZ              ) ;

	G4bool    IsPVPickable     ()                { return fPVPickable   ;}  
	void      SetPVPickability  ( G4bool on_off ) { fPVPickable = on_off ;}  
	G4double  SetPVTransparency ()  ; 
	G4double  GetPVTransparency () { return fPVTransparency ; } 

// data 
private:

	G4VRML2& fSystem;	// Graphics system for this scene.

	G4bool   fPVPickable     ;
	G4double fPVTransparency ;

	static G4int fSceneIdCount;
	static G4int fSceneCount;    // No. of existing scenes.

public: 
	G4FRClient fDest ;

};

#endif //G4VRML2_SCENE_HH
#endif //G4VIS_BUILD_VRML_DRIVER
