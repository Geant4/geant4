// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML2FileSceneHandler.hh,v 1.5 1999-12-16 17:25:14 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML2FileSceneHandler.hh
// Satoshi Tanaka & Yasuhide Sawada

#ifdef  G4VIS_BUILD_VRMLFILE_DRIVER

#ifndef G4VRML2FILE_SCENE_HANDLER_HH
#define G4VRML2FILE_SCENE_HANDLER_HH

#include "g4std/fstream"

#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4VSceneHandler.hh"


class G4VRML2File;
class G4VisAttributes;

class G4VRML2FileSceneHandler: public G4VSceneHandler {


// methods (public) 
public:
	G4VRML2FileSceneHandler(G4VRML2File& system, const G4String& name = "");
	virtual ~G4VRML2FileSceneHandler();
	void AddThis(const G4Box&);
	void AddThis(const G4Cons&);
	void AddThis(const G4Tubs&);
	void AddThis(const G4Trd&);
	void AddThis(const G4Trap&);
	void AddThis(const G4Sphere&);
        void AddThis(const G4Para&);
	void AddThis(const G4Torus&);
        void AddThis ( const G4Polycone& polycone ) {
          G4VSceneHandler::AddThis (polycone);
        }
        void AddThis ( const G4Polyhedra& polyhedra) {
          G4VSceneHandler::AddThis (polyhedra);
        }
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

	void connectPort();
	void closePort();

// methods (private) 
private:

	void      SendMaterialNode          ( const G4VisAttributes*  pAV ); 
	void      SendMaterialNode          ();

	void      SendLineColor             ( const G4VisAttributes*  pAV ); 
	void      SendMarkerColor           ( const G4VMarker&  mark ) ;
	void      SendMarkerWorldPosition   ( const G4VMarker&  mark ) ;

	G4double  GetMarkerHalfSize         ( const G4VMarker&  mark ) ;
	void      GetMarkerWorldPosition    (	const G4VMarker&  mark , 
						double* pX             ,
						double* pY             ,
						double* pZ              ) ;

	G4bool    isConnected      () { return fFlagDestOpen ; }

	G4bool    IsPVPickable     ()                { return fPVPickable   ;}  
	void      SetPVPickability ( G4bool on_off ) { fPVPickable = on_off ;}  
	G4double  SetPVTransparency ()  ; 
	G4double  GetPVTransparency () { return fPVTransparency ; } 

// data 
private:
	char fVRMLFileDestDir[256] ; 
	char fVRMLFileName[256]    ; 

	G4VRML2File& fSystem;	// Graphics system for this scene.
	G4bool       fFlagDestOpen ;

	G4int        fMaxFileNum   ;	

	G4bool       fPVPickable ;
	G4double     fPVTransparency ;

	static G4int fSceneIdCount;
	static G4int fSceneCount;    // No. of existing scenes.

public: 
	G4std::ofstream     fDest ;

};

#endif //G4VRML2File_SCENE_HH
#endif //G4VIS_BUILD_VRMLFILE_DRIVER
