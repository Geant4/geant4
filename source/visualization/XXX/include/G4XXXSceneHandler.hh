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
// $Id: G4XXXSceneHandler.hh,v 1.7 2001-11-21 16:47:32 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A base class for a scene handler to dump geometry hierarchy.

#ifndef G4XXXSCENEHANDLER_HH
#define G4XXXSCENEHANDLER_HH

#define G4XXXDEBUG  // Comment this out to suppress debug code.

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

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4ModelingParameters;

class G4XXXSceneHandler: public G4VSceneHandler {

public:
  G4XXXSceneHandler(G4VGraphicsSystem& system,
		      const G4String& name);
  virtual ~G4XXXSceneHandler();

  ////////////////////////////////////////////////////////////////
  // No need to implement these, but if you do...
  void AddThis(const G4Box&);
  void AddThis(const G4Cons&);
  void AddThis(const G4Tubs&);
  void AddThis(const G4Trd&);
  void AddThis(const G4Trap&);
  void AddThis(const G4Sphere&);
  void AddThis(const G4Para&);
  void AddThis(const G4Torus&);
  void AddThis(const G4Polycone&);
  void AddThis(const G4Polyhedra&);
  void AddThis(const G4VSolid&);
  // void PreAddThis(const G4Transform3D& objectTransformation,
  //                 const G4VisAttributes&);
  // void PostAddThis();

  // Allows G4PhysicalVolumeModel to provide extra information...
  void EstablishSpecials(G4PhysicalVolumeModel&);

  ////////////////////////////////////////////////////////////////
  // Required implementation of pure virtual functions...

  void AddPrimitive(const G4Polyline&);
  void AddPrimitive(const G4Text&);
  void AddPrimitive(const G4Circle&);
  void AddPrimitive(const G4Square&);
  void AddPrimitive(const G4Polyhedron&);
  void AddPrimitive(const G4NURBS&);

  ////////////////////////////////////////////////////////////////
  // Further optional AddPrimtive methods.  Explicitly invoke base
  // class methods if not otherwise defined to avoid warnings about
  // hiding of base class methods.
  void AddPrimitive(const G4Polymarker& polymarker) {
    G4VSceneHandler::AddPrimitive (polymarker);
  }
  void AddPrimitive(const G4Scale& scale) {
    G4VSceneHandler::AddPrimitive (scale);
  }

  ////////////////////////////////////////////////////////////////
  // Further optional AddPrimtive methods...

  // void AddPrimitive(const G4Scale&);
  // void AddPrimitive(const G4Polymarker&);

  ////////////////////////////////////////////////////////////////
  // Further optional virtual functions...

  // void BeginPrimitives(const G4Transform3D& objectTransformation);
  // void EndPrimitives();

  // void BeginModeling();
  // void EndModeling();

  ////////////////////////////////////////////////////////////////
  // Required...

  static G4int GetSceneCount() {return fSceneCount;}

protected:
  static G4int         fSceneIdCount;  // Counter for XXX scene handlers.
  static G4int         fSceneCount;    // No. of extanct scene handlers.
  G4int                fCurrentDepth;  // Current depth of geom. hierarchy.
  G4VPhysicalVolume*   fpCurrentPV;    // Current physical volume.
  G4LogicalVolume*     fpCurrentLV;    // Current logical volume.

private:
#ifdef G4XXXDEBUG
  void PrintThings();
#endif
};

#endif
