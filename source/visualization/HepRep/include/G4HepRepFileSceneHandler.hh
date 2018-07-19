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
// $Id: G4HepRepFileSceneHandler.hh 99152 2016-09-07 08:04:30Z gcosmo $
//
//
// Joseph Perl  27th January 2002
// A base class for a scene handler to export geometry and trajectories
// to the HepRep xml file format.

#ifndef G4HepRepFileSCENEHANDLER_HH
#define G4HepRepFileSCENEHANDLER_HH

//#define G4HEPREPFILEDEBUG  // Comment this out to suppress debug code.

#include "G4VSceneHandler.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Sphere.hh"
#include "G4Para.hh"
#include "G4Torus.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"

// HepRep
#include "G4HepRepFileXMLWriter.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4ModelingParameters;
class G4VisTrajContext;

class G4HepRepFileSceneHandler: public G4VSceneHandler {

public:
  G4HepRepFileSceneHandler(G4VGraphicsSystem& system,
		      const G4String& name);
  virtual ~G4HepRepFileSceneHandler();

  ////////////////////////////////////////////////////////////////
  // No need to implement these, but if you do...
  void AddSolid(const G4Box&);
  void AddSolid(const G4Cons&);
  void AddSolid(const G4Tubs&);
  void AddSolid(const G4Trd&);
  void AddSolid(const G4Trap&);
  void AddSolid(const G4Sphere&);
  void AddSolid(const G4Para&);
  void AddSolid(const G4Torus&);
  void AddSolid(const G4Polycone&);
  void AddSolid(const G4Polyhedra&);
  void AddSolid(const G4Orb&);
  void AddSolid(const G4Ellipsoid&);
  void AddSolid(const G4VSolid&);
  void AddCompound (const G4VTrajectory&);
  void InitTrajectory();
  void AddCompound (const G4VHit&);
  void InitHit();
  void AddCompound (const G4THitsMap<G4double>& hits) {
    G4VSceneHandler::AddCompound(hits);
  }
  void AddCompound (const G4THitsMap<G4StatDouble>& hits) {
    G4VSceneHandler::AddCompound(hits);
  }
  void AddCompound (const G4VDigi& digi) {
    G4VSceneHandler::AddCompound(digi);
  }
  // void PreAddSolid(const G4Transform3D& objectTransformation,
  //                 const G4VisAttributes&);
  // void PostAddSolid();

  ////////////////////////////////////////////////////////////////
  // Required implementation of pure virtual functions...

  void AddPrimitive(const G4Polyline&);
  void AddPrimitive(const G4Text&);
  void AddPrimitive(const G4Circle&);
  void AddPrimitive(const G4Square&);
  void AddPrimitive(const G4Polyhedron&);

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // Further optional AddPrimtive methods.  Explicitly invoke base
  // class methods if not otherwise defined to avoid warnings about
  // hiding of base class methods.
  void AddPrimitive(const G4Polymarker&);
  void AddPrimitive(const G4Scale& scale) {
    G4VSceneHandler::AddPrimitive (scale);
  }

  ////////////////////////////////////////////////////////////////
  // Further optional virtual functions...

  // void BeginPrimitives(const G4Transform3D& objectTransformation);
  // void EndPrimitives();

  void BeginModeling();
  void EndModeling();
  
  void BeginPrimitives2D(const G4Transform3D& objectTransformation);
  void EndPrimitives2D();

  //////////////////////////////////////////////////////////////
  // Administration functions.

  //void ClearStore ();
  void ClearTransientStore ();

  ////////////////////////////////////////////////////////////////
  // Required...

  G4HepRepFileXMLWriter *GetHepRepXMLWriter();

protected:
  static G4int         fSceneIdCount;  // Counter for HepRep scene handlers.

private:
  G4HepRepFileXMLWriter *hepRepXMLWriter;
  void AddHepRepInstance(const char* primName,
			 const G4Visible visible);
  void CheckFileOpen();
  int fileCounter;
  G4bool haveVisible;
  G4bool inPrimitives2D;
  G4bool warnedAbout3DText;
  G4bool warnedAbout2DMarkers;
  G4bool drawingTraj;
  G4bool doneInitTraj;
  G4bool drawingHit;
  G4bool doneInitHit;
  
  const G4VisTrajContext* trajContext;
  
  std::vector<G4AttValue>* trajAttValues;
  std::map<G4String,G4AttDef>* trajAttDefs;
  std::vector<G4AttValue>* hitAttValues;
  std::map<G4String,G4AttDef>* hitAttDefs;
  
#ifdef G4HEPREPFILEDEBUG
  void PrintThings();
#endif

};

#endif
