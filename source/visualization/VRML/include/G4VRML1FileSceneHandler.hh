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
//
//
// G4VRML1FileSceneHandler.hh
// Satoshi Tanaka & Yasuhide Sawada

#ifndef G4VRML1FILE_SCENE_HANDLER_HH
#define G4VRML1FILE_SCENE_HANDLER_HH

#include <fstream>

#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4VSceneHandler.hh"


class G4VRML1File;
class G4VisAttributes;

class G4VRML1FileSceneHandler: public G4VSceneHandler {

  friend class G4VRML1FileViewer;

// methods (public) 
public:
	G4VRML1FileSceneHandler(G4VRML1File& system, const G4String& name = "");
	virtual ~G4VRML1FileSceneHandler();
	void AddSolid(const G4Box&);
	void AddSolid(const G4Cons&);
	void AddSolid(const G4Tubs&);
	void AddSolid(const G4Trd&);
	void AddSolid(const G4Trap&);
	void AddSolid(const G4Sphere&);
        void AddSolid(const G4Para&);
	void AddSolid(const G4Torus&);
        void AddSolid ( const G4Polycone& polycone ) {
          G4VSceneHandler::AddSolid (polycone);
        }
        void AddSolid ( const G4Polyhedra& polyhedra) {
          G4VSceneHandler::AddSolid (polyhedra);
        }
        void AddSolid ( const G4Orb& orb ) {
          G4VSceneHandler::AddSolid (orb);
        }
        void AddSolid ( const G4Ellipsoid& ellipsoid) {
          G4VSceneHandler::AddSolid (ellipsoid);
        }
        void AddSolid ( const G4TessellatedSolid& tess) {
          G4VSceneHandler::AddSolid (tess);
        }
        void AddSolid(const G4VSolid&);
        void AddCompound ( const G4VTrajectory& traj) {
          G4VSceneHandler::AddCompound(traj);
        }
        void AddCompound ( const G4VHit& hit) {
          G4VSceneHandler::AddCompound(hit);
        }
        void AddCompound ( const G4VDigi& digi) {
          G4VSceneHandler::AddCompound(digi);
        }
        void AddCompound ( const G4THitsMap<G4double> & hits) {
	  G4VSceneHandler::AddCompound(hits);
	}
        void AddCompound ( const G4THitsMap<G4StatDouble> & hits) {
	  G4VSceneHandler::AddCompound(hits);
	}

	void BeginPrimitives(const G4Transform3D& objectTransformation);
	void EndPrimitives();

	void AddPrimitive(const G4Polyline&);
	void AddPrimitive(const G4Polyhedron&);
	void AddPrimitive(const G4Text&); 
	void AddPrimitive(const G4Circle&);
	void AddPrimitive(const G4Square&);
	void AddPrimitive (const G4Polymarker& polymarker)
		{ G4VSceneHandler::AddPrimitive (polymarker); }
        void AddPrimitive (const G4Scale& scale) 
                { G4VSceneHandler::AddPrimitive (scale); }

	void ClearTransientStore();  // Used for triggering detector re-drawing.

	void BeginModeling();
	void EndModeling();

	void VRMLBeginModeling();
	void VRMLEndModeling();

	void connectPort();
	void closePort();

// methods (private) 
private:

	void SendMaterialNode( const G4VisAttributes*  pAV ); // NEW NEW NEW
	void SendMaterialNode();
	void SendMatrixTransformNode(const G4Transform3D &);
	void SendCubeNode(G4double, G4double, G4double);
	void SendCylinderNode(G4double, G4double);
	void SendSphereNode(G4double);

	void      SendMarkerColor         ( const G4VMarker&  mark ) ;
	void      SendMarkerWorldPosition ( const G4VMarker&  mark ) ;
	G4double  GetMarkerHalfSize       ( const G4VMarker&  mark ) ;

	G4bool    isConnected      () { return fFlagDestOpen ; }

// data 
private:
	char fVRMLFileDestDir[256] ; 
	char fVRMLFileName[256]    ; 

	G4String fCurrentDEF;

	G4VRML1File& fSystem;	// Graphics system for this scene.
	std::ofstream     fDest ;
	G4bool       fFlagDestOpen ;

	G4int        fMaxFileNum   ;	

	static G4int fSceneIdCount;

};

#endif //G4VRML1FILE_SCENE_HANDLER_HH

