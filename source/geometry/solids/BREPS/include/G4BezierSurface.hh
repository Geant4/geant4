// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BezierSurface.hh,v 1.5 2000-11-08 14:22:00 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4BezierSurface
//
// Class description:
// 
// Definition of a generic Bezier surface.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __BEZIERSURFACE_H
#define __BEZIERSURFACE_H

#include "G4Ray.hh"
#include "G4ControlPoints.hh"
#include "G4SurfaceList.hh"
#include "G4PointRat.hh"
#include "G4OsloMatrix.hh"
#include "G4KnotVector.hh"     

class G4ProjectedSurface;

class G4BezierSurface : public G4Surface
{
  friend class G4BSplineSurface;
  friend class G4ProjectedSurface; 
   
public:  // with description
  
  G4BezierSurface();
  virtual ~G4BezierSurface();
    // Constructor & destructor

  inline G4Point3D AveragePoint() const;
  inline void SetAveragePoint(G4Point3D p);

  inline G4double UAverage() const;
  inline G4double VAverage() const;

  inline void Dir(G4int d);
  inline void ChangeDir();

  inline G4double SMin() const;
  inline G4double SMax() const;

  inline G4int GetOrder(G4int direction) const;
  inline void PutOrder(G4int direction, G4int value);

  inline G4double GetU() const;
  inline G4double GetV() const;

  void CalcBBox();

  G4int BIntersect(G4SurfaceList&);

  G4int ClipBothDirs();
  void ClipSurface();

  virtual G4Vector3D SurfaceNormal(const G4Point3D& Pt) const;

  friend void CopySurface(G4BezierSurface& bez);

protected:  // without description

  G4SurfaceList* bezier_list;

  static G4int Clips;
  static G4int Splits;    
  static G4double Tolerance;
    // Test variables 

private:

  void CalcAverage();
  void CalcDistance(const G4Point3D&);
  void SetValues();
  inline void LocalizeClipValues();

  void SplitNURBSurface();
  void GetClippedRegionFromSurface();
  void RefineSurface();
  void CalcOsloMatrix();
  void MapSurface(G4Surface*);

  inline G4double Findzero(G4double x0, G4double x1,
                           G4double y0, G4double y1) const;
  inline G4int Sign(G4double a) const;
    // For ClipSurface...

  inline G4int Amax(G4int i, G4int j) const;
  inline G4int Amin(G4int i, G4int j) const;
  inline G4int AhIndex(G4int j, G4int t, G4int iorder) const;
    // For calc_G4OsloMatrix...

private:

  G4BezierSurface(const G4BezierSurface&);
  G4BezierSurface& operator=(const G4BezierSurface &);
    // Private copy constructor and assignment operator.

private:

  G4int           order[2];
  G4double        smin;
  G4double        smax;
  G4Point3D       line;
  G4double        average_u;
  G4double        average_v;
  G4Point3D       average_pt;
  G4int           dir;    
  G4KnotVector    *u_knots;
  G4KnotVector    *v_knots;
  G4ControlPoints *ctl_points;
  
  G4KnotVector     *new_knots;
  G4int            ord;
  G4OsloMatrix     * oslo_m;
  G4int            lower,upper;
  G4double         u[2];
  G4double         v[2];
  G4double         u_min;
  G4double         u_max;
  G4double         v_min;
  G4double         v_max;
  G4ControlPoints* old_points;  

};

#include "G4BezierSurface.icc"

#endif
