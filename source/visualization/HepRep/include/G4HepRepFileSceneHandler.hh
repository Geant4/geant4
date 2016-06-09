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
// $Id: G4HepRepFileSceneHandler.hh,v 1.22 2005/06/02 17:43:46 allison Exp $
// GEANT4 tag $Name: geant4-07-01 $
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
  void AddSolid(const G4VSolid&);
  void AddCompound (const G4VTrajectory&);
  void AddCompound (const G4VHit&);
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
  void AddPrimitive(const G4NURBS&);

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
  char fileDir[256]; 
  char fileName[256];
  G4bool fileOverwrite;
  G4bool cullInvisibleObjects;
  G4bool haveVisible;
  G4bool drawingTraj;
  G4bool drawingHit;

#ifdef G4HEPREPFILEDEBUG
  void PrintThings();
#endif

};

#endif
