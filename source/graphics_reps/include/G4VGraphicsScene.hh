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
// $Id: G4VGraphicsScene.hh 102801 2017-02-22 15:17:53Z gcosmo $
// John Allison  19th July 1996
//
// Class Description:
// Abstract interface class for a graphics scene handler.
// It is a minimal scene handler for the GEANT4 kernel.
// See G4VSceneHandler for a fuller description.  G4VSceneHandler is
// the full abstract interface to graphics systems.

#ifndef G4VGRAPHICSSCENE_HH
#define G4VGRAPHICSSCENE_HH

#include "globals.hh"
#include "G4Transform3D.hh"

class G4VisAttributes;
class G4VSolid;
class G4Box;
class G4Cons;
class G4Orb;
class G4Para;
class G4Torus;
class G4Trap;
class G4Trd;
class G4Tubs;
class G4Sphere;

class G4Ellipsoid;
class G4Polycone;
class G4Polyhedra;

class G4PhysicalVolumeModel;
class G4VTrajectory;
class G4VHit;
class G4VDigi;
template <typename T> class G4THitsMap;
class G4Polyline;
class G4Scale;
class G4Text;
class G4Circle;
class G4Square;
class G4Polymarker;
class G4Polyhedron;
class G4VisExtent;
class G4StatDouble;

class G4VGraphicsScene {

public: // With description

  G4VGraphicsScene();
  virtual ~G4VGraphicsScene();

  ///////////////////////////////////////////////////////////////////
  // Methods for adding solids to the scene handler.  They
  // must always be called in the triplet PreAddSolid, AddSolid and
  // PostAddSolid.  The transformation and visualization attributes
  // must be set by the call to PreAddSolid.  A possible default
  // implementation is to request the solid to provide a G4Polyhedron
  // or similar primitive - see, for example, G4VSceneHandler in the
  // Visualization Category.

  virtual void PreAddSolid (const G4Transform3D& objectTransformation,
                            const G4VisAttributes& visAttribs) = 0;
  // objectTransformation is the transformation in the world
  // coordinate system of the object about to be added, and
  // visAttribs is its visualization attributes.

  virtual void PostAddSolid () = 0;

  // From geometry/solids/CSG
  virtual void AddSolid (const G4Box&)       = 0;
  virtual void AddSolid (const G4Cons&)      = 0;
  virtual void AddSolid (const G4Orb&)       = 0;
  virtual void AddSolid (const G4Para&)      = 0;
  virtual void AddSolid (const G4Sphere&)    = 0;
  virtual void AddSolid (const G4Torus&)     = 0;
  virtual void AddSolid (const G4Trap&)      = 0;
  virtual void AddSolid (const G4Trd&)       = 0;
  virtual void AddSolid (const G4Tubs&)      = 0;

  // From geometry/solids/specific
  virtual void AddSolid (const G4Ellipsoid&) = 0;
  virtual void AddSolid (const G4Polycone&)  = 0;
  virtual void AddSolid (const G4Polyhedra&) = 0;

  // For solids not above.
  virtual void AddSolid (const G4VSolid&)    = 0;

  ///////////////////////////////////////////////////////////////////
  // Methods for adding "compound" GEANT4 objects to the scene
  // handler.  These methods may either (a) invoke "user code" that
  // uses the "user interface", G4VVisManager (see, for example,
  // G4VSceneHandler in the Visualization Category, which for
  // trajectories uses G4VTrajectory::DrawTrajectory, via
  // G4TrajectoriesModel in the Modeling Category) or (b) invoke
  // AddPrimitives below (between calls to Begin/EndPrimitives) or (c)
  // use graphics-system-specific code or (d) any combination of the
  // above.

  virtual void AddCompound (const G4VTrajectory&)        = 0;
  virtual void AddCompound (const G4VHit&)               = 0;
  virtual void AddCompound (const G4VDigi&)              = 0;
  virtual void AddCompound (const G4THitsMap<G4double>&) = 0;
  virtual void AddCompound (const G4THitsMap<G4StatDouble>&) = 0;

  ///////////////////////////////////////////////////////////////////
  // Methods for adding graphics primitives to the scene handler.  A
  // sequence of calls to AddPrimitive must be sandwiched between
  // calls to BeginPrimitives and EndPrimitives.  A sequence is any
  // number of calls that have the same transformation.

  virtual void BeginPrimitives
  (const G4Transform3D& objectTransformation = G4Transform3D()) = 0;
  // objectTransformation is the transformation in the world
  // coordinate system of the object about to be added.

  virtual void EndPrimitives () = 0;

  virtual void BeginPrimitives2D
  (const G4Transform3D& objectTransformation = G4Transform3D()) = 0;

  virtual void EndPrimitives2D () = 0;
  // The x,y coordinates of the primitives passed to AddPrimitive are
  // intrepreted as screen coordinates, -1 < x,y < 1.  The
  // z-coordinate is ignored.

  virtual void AddPrimitive (const G4Polyline&)   = 0;
  virtual void AddPrimitive (const G4Scale&)      = 0;
  virtual void AddPrimitive (const G4Text&)       = 0;
  virtual void AddPrimitive (const G4Circle&)     = 0;
  virtual void AddPrimitive (const G4Square&)     = 0;
  virtual void AddPrimitive (const G4Polymarker&) = 0;
  virtual void AddPrimitive (const G4Polyhedron&) = 0;

  virtual const G4VisExtent& GetExtent() const;
  // The concrete class should overload this or a null extent will be returned.
  // See G4VScenHandler for example.

};

#endif
