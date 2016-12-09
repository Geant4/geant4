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
// $Id: G4GMocrenFileSceneHandler.hh 99152 2016-09-07 08:04:30Z gcosmo $
//
//
// Created:  Mar. 31, 2009  Akinori Kimura  
//
// Scene handler to export geometry and trajectories to a gMocren file.
//
#ifndef G4GMocrenFile_SCENE_HANDLER_HH
#define G4GMocrenFile_SCENE_HANDLER_HH

#include "globals.hh"

#include "G4VSceneHandler.hh"

#include <fstream>

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
  void AddSolid ( const G4Orb& orb ) {
    G4VSceneHandler::AddSolid (orb);
  }
  void AddSolid ( const G4Ellipsoid& ellipsoid) {
    G4VSceneHandler::AddSolid (ellipsoid);
  }
  void AddSolid ( const G4VSolid& solid  );
  void AddCompound ( const G4VTrajectory& traj);
  void AddCompound ( const G4VHit& hit);
  void AddCompound ( const G4VDigi& hit);
  void AddCompound ( const G4THitsMap<G4double> & hits);
  void AddCompound ( const G4THitsMap<G4StatDouble> & hits);


  void ClearTransientStore();  // Used for triggering detector re-drawing.

  //----- public methods inherent to this class
  void         GFBeginModeling () ;
  void         GFEndModeling   () ;
  G4bool       GFIsInModeling  () { return kFlagInModeling ; }

  G4bool IsSavingGdd   ( void ) { return kFlagSaving_g4_gdd ;	}
  void	BeginSavingGdd( void ); 
  void	EndSavingGdd  ( void ) ;
  void	SetGddFileName() ;

  G4GMocrenFile&  GetSystem   () { return kSystem   ; }
  const char*  GetGddFileName () { return kGddFileName ; }


private:

  //----- initialize all parameters
  void InitializeParameters();

  //----- Utilities etc.
  G4bool IsVisible();

  //
  void AddDetector(const G4VSolid & solid);
  void ExtractDetector();

  void GetNestedVolumeIndex(G4int, G4int[3]);

private:
  G4GMocrenFile&	kSystem;     // Graphics system for this scene.
  G4GMocrenMessenger & kMessenger;
  G4GMocrenIO * kgMocrenIO;

  std::map<G4int, float> kModality;
  G4int kModalitySize[3];
  //std::map<G4ThreeVector, float> kModalityDensities; // key: position, val: density
  G4bool kbSetModalityVoxelSize;
  G4bool kbModelingTrajectory;

  static G4int	kSceneIdCount;
  //std::vector<float *> fTrajectories;
  //std::vector<unsigned char *> fTrajectoryColors;
  G4Transform3D kVolumeTrans3D;

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
  std::vector<Detector> kDetectors;
  G4ThreeVector kVolumeSize;
  G4ThreeVector kVoxelDimension;
  std::vector<G4String> kNestedVolumeNames;
  G4int kNestedVolumeDimension[3];
  G4int kNestedVolumeDirAxis[3];

  class Index3D {
  public:
    G4int x, y, z;

    Index3D();
    Index3D(const Index3D & _index3D);
    Index3D(G4int _x, G4int _y, G4int _z);
    ~Index3D(){;}
    G4bool operator < (const Index3D & _right) const;
    G4bool operator == (const Index3D & _right) const;
  private:
    // Private assigment operator -
    // assignment not allowed.  Keeps Coverity happy.
    // Index3D& operator = (const Index3D&);
  };

  std::map<Index3D, float> kNestedModality;
  //std::map<Index3D, G4double> * fTempNestedHits;
  std::map<G4String, std::map<Index3D, G4double> > kNestedHitsList;
  //std::map<G4String, G4String> kNestedHitsUnit;

  std::ofstream	kGddDest;  // defined here
  G4bool kFlagInModeling;	
		// true:  GF_BEGIN_MODELING has sent to gMocrenFile, and
		//        GF_END_MODELING   has not sent yet.
		// false:  otherwise
		// 
		// kFlagInModeling is set to "true"
		// in GFBeginModeling(), and to "false" 
		// in GFEndModeling().

  G4bool kFlagSaving_g4_gdd ;	

  G4int kFlagParameterization; // 0: G4VNestedParameterisation based geometry
                               // 1: G4PhantomParameterisation
                               // 2: interactive scorer
  G4bool kFlagProcessedInteractiveScorer;

  char kGddDestDir[256]; 
  char kGddFileName[256];
  G4int	kMaxFileNum;

};

#endif
