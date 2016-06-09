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
// $Id: G4BoundingSphereScene.hh,v 1.20 2010-05-30 11:23:25 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  7th June 1997
// An artificial scene to reuse G4VScene code to calculate a bounding sphere.

#ifndef G4BOUNDINGSPHERESCENE_HH
#define G4BOUNDINGSPHERESCENE_HH

#include "G4VGraphicsScene.hh"
#include "G4VisExtent.hh"
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

class G4VModel;

class G4BoundingSphereScene: public G4VGraphicsScene {

public:
  G4BoundingSphereScene (G4VModel* pModel = 0);
  virtual ~G4BoundingSphereScene ();
  void PreAddSolid (const G4Transform3D& objectTransformation,
		    const G4VisAttributes&);
  void PostAddSolid () {}
  void AddSolid (const G4Box& s) {Accrue (s);}
  void AddSolid (const G4Cons& s) {Accrue (s);}
  void AddSolid (const G4Tubs& s) {Accrue (s);}
  void AddSolid (const G4Trd& s) {Accrue (s);}
  void AddSolid (const G4Trap& s) {Accrue (s);}
  void AddSolid (const G4Sphere& s) {Accrue (s);}
  void AddSolid (const G4Para& s) {Accrue (s);}
  void AddSolid (const G4Torus& s) {Accrue (s);}
  void AddSolid (const G4Polycone& s) {Accrue (s);}
  void AddSolid (const G4Polyhedra& s) {Accrue (s);}
  void AddSolid (const G4VSolid& s) {Accrue (s);}
  void AddCompound (const G4VTrajectory&) {}
  void AddCompound (const G4VHit&) {}
  void AddCompound (const G4VDigi&) {}
  void AddCompound (const G4THitsMap<G4double>&) {}
  G4VisExtent GetBoundingSphereExtent ();
  const G4Point3D& GetCentre() const {return fCentre;}
  G4double GetRadius() const {return fRadius;}

  void SetCentre(const G4Point3D& centre) {fCentre = centre;}

  ////////////////////////////////////////////////////////////////
  // The following 2 functions can be used by any code which wishes to
  // accrue a bounding sphere.  Just instantiate a
  // G4BoundingSphereScene and use AccrueBoundingSphere.

  void ResetBoundingSphere ();
  void AccrueBoundingSphere (const G4Point3D& centre,
			     G4double radius);

  ////////////////////////////////////////////////////////////////
  // Functions not used by required by the abstract interface.

  virtual void BeginPrimitives (const G4Transform3D&) {}
  virtual void EndPrimitives () {}
  virtual void BeginPrimitives2D (const G4Transform3D&) {}
  virtual void EndPrimitives2D () {}
  virtual void AddPrimitive (const G4Polyline&)   {}
  virtual void AddPrimitive (const G4Scale&)      {}
  virtual void AddPrimitive (const G4Text&)       {}
  virtual void AddPrimitive (const G4Circle&)     {}
  virtual void AddPrimitive (const G4Square&)     {}
  virtual void AddPrimitive (const G4Polymarker&) {}
  virtual void AddPrimitive (const G4Polyhedron&) {}
  virtual void AddPrimitive (const G4NURBS&)      {}

private:
  void Accrue (const G4VSolid& solid);
  G4VModel* fpModel;  // Instantiating code may optionally set this.
  G4Point3D fCentre;
  G4double fRadius;
  const G4Transform3D* fpObjectTransformation;
};

#endif
