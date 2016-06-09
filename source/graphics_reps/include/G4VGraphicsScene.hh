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
// $Id: G4VGraphicsScene.hh,v 1.5 2005/11/10 15:39:16 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $
// John Allison  19th July 1996
//
// Class Description:
// Abstract interface class for a graphics scene handler.
// It is a minimal scene handler for the GEANT4 kernel.
// See G4VSceneHandler for a fuller description.  G4VSceneHandler is
// the full abstract interface to graphics systems.

#ifndef G4VGRAPHICSSCENE_HH
#define G4VGRAPHICSSCENE_HH

class G4VisAttributes;
class G4VSolid;
class G4Box;
class G4Cons;
class G4Tubs;
class G4Trd;
class G4Trap;
class G4Sphere;
class G4Ellipsoid;
class G4Para;
class G4Torus;
class G4PhysicalVolumeModel;
class G4Polycone;
class G4Polyhedra;
class G4VTrajectory;
class G4VHit;
class G4Polyline;
class G4Scale;
class G4Text;
class G4Circle;
class G4Square;
class G4Polymarker;
class G4Polyhedron;
class G4NURBS;

#include "G4Transform3D.hh"

class G4VGraphicsScene {

public: // With description

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

  virtual void AddSolid (const G4Box&)       = 0;
  virtual void AddSolid (const G4Cons&)      = 0;
  virtual void AddSolid (const G4Tubs&)      = 0;
  virtual void AddSolid (const G4Trd&)       = 0;
  virtual void AddSolid (const G4Trap&)      = 0;
  virtual void AddSolid (const G4Sphere&)    = 0;
  virtual void AddSolid (const G4Para&)      = 0;
  virtual void AddSolid (const G4Torus&)     = 0;
  virtual void AddSolid (const G4Polycone&)  = 0;
  virtual void AddSolid (const G4Polyhedra&) = 0;
  virtual void AddSolid (const G4VSolid&)    = 0;  // For solids not above.

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

  virtual void AddCompound (const G4VTrajectory&) = 0;
  virtual void AddCompound (const G4VHit&)        = 0;

  ///////////////////////////////////////////////////////////////////
  // Methods for adding graphics primitives to the scene handler.  A
  // sequence of calls to AddPrimitive must be sandwiched between
  // calls to BeginPrimitives and EndPrimitives.  A sequence is any
  // number of calls the have the same transformation.

  virtual void BeginPrimitives
  (const G4Transform3D& objectTransformation = G4Transform3D()) = 0;
  // objectTransformation is the transformation in the world
  // coordinate system of the object about to be added.

  virtual void EndPrimitives () = 0;

  virtual void AddPrimitive (const G4Polyline&)   = 0;
  virtual void AddPrimitive (const G4Scale&)      = 0;
  virtual void AddPrimitive (const G4Text&)       = 0;
  virtual void AddPrimitive (const G4Circle&)     = 0;
  virtual void AddPrimitive (const G4Square&)     = 0;
  virtual void AddPrimitive (const G4Polymarker&) = 0;
  virtual void AddPrimitive (const G4Polyhedron&) = 0;
  virtual void AddPrimitive (const G4NURBS&)      = 0;

  ///////////////////////////////////////////////////////////////////
  // Special functions for particular models.

  virtual void EstablishSpecials (G4PhysicalVolumeModel&) {}
  // Used to establish any special relationships between scene and
  // this particular type of model - non-pure, i.e., no requirement to
  // implement.  See G4PhysicalVolumeModel.hh for details.

  virtual void DecommissionSpecials (G4PhysicalVolumeModel&) {}
  // Used to reverse the effect of EstablishSpecials, if required.

};

#endif
