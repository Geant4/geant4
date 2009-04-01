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
// $Id: G4GMocrenFileSceneHandler.hh,v 1.1 2009-04-01 13:16:11 akimura Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Akinori Kimura    March 31, 2009
//
// Scene handler to export geometry and trajectories to a gMocren file.
//
#ifndef G4GMocrenFile_SCENE_HANDLER_HH
#define G4GMocrenFile_SCENE_HANDLER_HH

#include "globals.hh"

#include "G4VSceneHandler.hh"

#include "G4FRofstream.hh"
#include "G4FRConst.hh"

#include "G4THitsMap.hh"

class G4VisAttributes ;
class G4GMocrenFile;
class G4GMocrenMessenger;
class G4GMocrenIO;
class G4VSolid;
class G4Polyhedron;
class G4Colour;

	//-----
class G4GMocrenFileSceneHandler: public G4VSceneHandler {

  friend class G4GMocrenFileViewer;

public:

	//----- constructor and destructor
  G4GMocrenFileSceneHandler (G4GMocrenFile& system,
			     G4GMocrenMessenger & messenger,
			     const G4String& name = "");
  virtual ~G4GMocrenFileSceneHandler ();

	//----- overriding base class methods
  void AddPrimitive (const G4Polyline& line);
  void AddPrimitive (const G4Polyhedron& p);
  void AddPrimitive (const G4NURBS& nurb);
  void AddPrimitive (const G4Text&);
  void AddPrimitive (const G4Circle&);
  void AddPrimitive (const G4Square&);

	//----- explicitly invoke base class methods to avoid warnings about
        //----- hiding of base class methods.
  void AddPrimitive (const G4Polymarker& polymarker) 
       { G4VSceneHandler::AddPrimitive (polymarker); }
  void AddPrimitive (const G4Scale& scale) 
       { G4VSceneHandler::AddPrimitive (scale); }

  virtual void BeginModeling () { G4VSceneHandler::BeginModeling ();} 
  virtual void EndModeling   () { G4VSceneHandler::EndModeling   ();}

  virtual void BeginPrimitives (const G4Transform3D& objectTransformation);
  virtual void EndPrimitives ();

  void AddSolid ( const G4Box&    box    );
  void AddSolid ( const G4Cons&   cons   );
  void AddSolid ( const G4Tubs&   tubs   );
  void AddSolid ( const G4Trd&    trd    );
  void AddSolid ( const G4Trap&   trap   );
  void AddSolid ( const G4Sphere& sphere );
  void AddSolid ( const G4Para&   para   );
  void AddSolid ( const G4Torus&  torus  );
  void AddSolid ( const G4Polycone& polycone ) {
    G4VSceneHandler::AddSolid (polycone);
  }
  void AddSolid ( const G4Polyhedra& polyhedra) {
    G4VSceneHandler::AddSolid (polyhedra);
  }
  void AddSolid ( const G4VSolid& solid  );
  void AddCompound ( const G4VTrajectory& traj);
  void AddCompound ( const G4VHit& hit);
  void AddCompound ( const G4THitsMap<G4double> & hits);


  void ClearTransientStore();  // Used for triggering detector re-drawing.

  //----- public methods inherent to this class
  void         FRBeginModeling () ;
  void         FREndModeling   () ;
  G4bool       FRIsInModeling  () { return FRflag_in_modeling ; }

  G4bool IsSavingGdd   ( void ) { return flag_saving_g4_gdd ;	}
  void	BeginSavingGdd( void ); 
  void	EndSavingGdd  ( void ) ;
  void	SetGddFileName() ;

  G4GMocrenFile&  GetSystem   () { return fSystem   ; }
  const char*  GetGddFileName () { return fGddFileName ; }


private:

  //----- Utilities etc (common to DAWN and GMocrenFile drivers )
  G4bool    IsVisible     ( void ) ;

  //
  void AddDetector(const G4VSolid & solid);
  void ExtractDetector();

private:
  G4GMocrenFile&	fSystem;     // Graphics system for this scene.
  G4GMocrenMessenger & fMessenger;
  G4GMocrenIO * fgMocrenIO;

  std::map<G4int, float> fModality;
  short fCTMapMinMax[2];
  std::vector<float> fCTMap;
  G4int fModalitySize[3];
  //std::map<G4ThreeVector, float> fModalityDensities; // key: position, val: density
  G4bool fbSetModalityVoxelSize;
  G4bool fbModelingTrajectory;

  static G4int	fSceneIdCount;
  //std::vector<float *> fTrajectories;
  //std::vector<unsigned char *> fTrajectoryColors;
  G4Transform3D fVolumeTrans3D;

  class Detector {
  public:
    G4String name;
    G4Polyhedron * polyhedron;
    G4Transform3D transform3D;
    unsigned char color[3];
    Detector();
    ~Detector();
    void clear();
  };
  std::vector<Detector> fDetectors;
  G4ThreeVector fVolMinPosition, fVolMaxPosition;
  G4ThreeVector fVolumeSize;
  G4ThreeVector fVoxelDimension;
  std::vector<G4String> fNestedVolumeNames;

  class Index3D{
  public:
    G4int x, y, z;
    Index3D();
    Index3D(G4int _x, G4int _y, G4int _z);
    ~Index3D(){;}
    G4bool operator < (const Index3D & _right) const;
    G4bool operator == (const Index3D & _right) const;
  };
  std::map<Index3D, float> fNestedModality;
  std::map<Index3D, G4double> * fTempNestedHits;
  std::map<G4String, std::map<Index3D, G4double> > fNestedHitsList;
  //std::map<G4String, G4String> fNestedHitsUnit;

  G4FRofstream	fGddDest    ;  // defined here
  G4bool	FRflag_in_modeling ;	
		// true:  FR_BEGIN_MODELING has sent to DAWN, and
		//        FR_END_MODELING   has not sent yet.
		// false:  otherwise
		// 
		// FRflag_in_modeling is set to "true"
		// in FRBeginModeling(), and to "false" 
		// in FREndModeling().

  G4bool	flag_saving_g4_gdd ;	

  const int	COMMAND_BUF_SIZE    ;

  char		fGddDestDir [256] ; 
  char          fGddFileName[256] ;
  G4int		fMaxFileNum           ;

  G4int         fPrec, fPrec2 ;
	
};

#endif
