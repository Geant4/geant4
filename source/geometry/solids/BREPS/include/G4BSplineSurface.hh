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
// $Id: G4BSplineSurface.hh,v 1.11 2006-06-29 18:38:11 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4BSplineSurface
//
// Class description:
// 
// Definition of a generic BSpline surface.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __BSPLINESURFACE_H
#define __BSPLINESURFACE_H

#include "G4Point3D.hh"
#include "G4PointRat.hh"
#include "G4Surface.hh"
#include "G4ProjectedSurface.hh"
#include "G4UVHit.hh"

class G4BSplineSurface : public G4Surface
{

 public:  // with description

  G4BSplineSurface();
  G4BSplineSurface(const char* nurbfilename, G4Ray& rayref);
  G4BSplineSurface(G4int u, G4int v, G4KnotVector& u_kv, G4KnotVector& v_kv, 
		   G4ControlPoints& cp);    
  virtual ~G4BSplineSurface();
    // Constructors & destructor.

  G4int Intersect(const G4Ray&);
  void CalcBBox();
    // Finds the bounds of the b-spline surface.
    // The bounding box is used for a preliminary check of intersection.
  
  inline G4double GetUHit() const;
  inline G4double GetVHit() const;
    	 
  G4double ClosestDistanceToPoint(const G4Point3D&);

  inline void Reset();

  inline G4int GetRows() const;
  inline G4int GetCols() const;
  inline G4Point3D GetControlPoint(G4int a, G4int b) const;
    // Accessors for control points.

public:
 
  inline G4int MyType() const;

private:

  G4BSplineSurface(const G4BSplineSurface&);
  G4BSplineSurface& operator=(const G4BSplineSurface&);
    // Private copy constructor and assignment operator.

  void FindIntersections(const G4Ray&);

  inline G4int GetOrder(G4int direction) const;
  inline void PutOrder(G4int direction, G4int value);

  void AddHit(G4double u, G4double v);
  void ProjectNURBSurfaceTo2D( const G4Plane& ,const G4Plane&,
			       register G4ProjectedSurface*);
    // Projects the nurb surface so that the z-axis = ray. 

  G4ProjectedSurface* CopyToProjectedSurface(const G4Ray&);
  G4Point3D  FinalIntersection();

  // L. Broglia
  // Because  G4BSplineSurface::Evaluate hides the virtual function 
  // G4Surface::Evaluate(const G4Ray&), I modified the function name
  // G4Point3D  Evaluate();  
  G4Point3D  BSEvaluate();

  G4PointRat& InternalEvalCrv(G4int i, G4ControlPoints *crv);
  
  G4Point3D   Evaluation(const G4Ray&);

  inline G4Vector3D  SurfaceNormal(const G4Point3D& Pt) const;
  
private:
  
  G4SurfaceList bezier_list;
  G4SurfaceList projected_list;
  short dir;
  G4int order[2];
  G4KnotVector *u_knots;
  G4KnotVector *v_knots;
  G4KnotVector *tmp_knots;
  G4ControlPoints *ctl_points;
  G4UVHit* Hit;
  G4UVHit* first_hit;
  G4int ord;
  G4int k_index;
  G4double param;
  G4int Rational;
}; 

#include "G4BSplineSurface.icc"

#endif
