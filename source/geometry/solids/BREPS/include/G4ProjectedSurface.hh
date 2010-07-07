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
// $Id: G4ProjectedSurface.hh,v 1.8 2010-07-07 14:45:31 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ProjectedSurface
//
// Class description:
// 
// Definition of a projected surface.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4ProjectedSurface_h
#define __G4ProjectedSurface_h 1

#include "G4BezierSurface.hh"

class G4ProjectedSurface : public G4Surface
{
  friend class G4BSplineSurface;
  friend void CopySurface(G4ProjectedSurface& proj);
  
 public:  // with description

  G4ProjectedSurface();
  virtual ~G4ProjectedSurface();
    // Default constructor and destructor.

  void CalcBBox();
    // Finds the bounds of the 2D-projected nurb, it
    // calculates the bounds for a bounding rectangle
    // to the surface. The bounding rectangle is used
    // for a preliminary check of intersection.
        
 public:  // without description

  inline G4Vector3D SurfaceNormal(const G4Point3D& Pt) const;
    // Returns normal to surface (G4Vector3D(0,0,0)).

 private:

  G4ProjectedSurface(const G4ProjectedSurface&);
  G4ProjectedSurface& operator=(const G4ProjectedSurface&);
    // Private copy constructor and assignment operator.

  void CopySurface();
    // Copies the projected surface into a bezier surface
    // and adds it to the List of bezier surfaces.

  void ConvertToBezier (G4SurfaceList& p, G4SurfaceList& b);
    // Converts surface into a Bezier surface to b.

  inline G4int GetOrder(G4int direction) const;
  inline void PutOrder(G4int direction, G4int value);
  
  void SplitNURBSurface();
    // Divides the surface in two parts. Uses the oslo-algorithm to calculate
    // the new knot-vectors and control-points for the subsurfaces.

  G4int  CheckBezier();
    // Checks if the surface is a Bezier surface by verifying
    // if internal knots exist. If no internal knots exist the quantity
    // of knots is 2*order of the surface. Returns 1 if the surface 
    // is a Bezier.

  void CalcOsloMatrix();
    // Calculates the oslo-matrix, which is used in mapping the new
    // knot-vector and the control-point values.
    // This algorithm is described in the paper "Making the Oslo-algorithm
    // more efficient" in SIAM J.NUMER.ANAL. Vol.23, No. 3, June '86.

  void MapSurface(G4ProjectedSurface* srf);
    // Maps the new control-points into the new surface.
    // This algorithm is described in the paper "Making the Oslo-algorithm
    // more efficient" in SIAM J.NUMER.ANAL. Vol.23, No. 3, June '86.

  inline G4int Amax(G4int i, G4int j) const;
  inline G4int Amin(G4int i, G4int j) const;
  inline G4int AhIndex(G4int j, G4int t, G4int iorder) const;

 protected:

  static G4int Splits;
  G4ControlPoints *ctl_points;
    // Test variables

 private:

  short dir;
  G4KnotVector *u_knots;
  G4KnotVector *v_knots;
  
  G4SurfaceList* projected_list;
  G4SurfaceList* bezier_list;   
	
  G4int order[2];
  G4KnotVector *new_knots;
  G4int ord;
  G4int lower,upper;
  
  G4OsloMatrix* oslo_m;
  G4Point3D vmin;
  G4Point3D vmax;
};	

#include "G4ProjectedSurface.icc"

#endif
