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
// $Id: G4BSplineCurve.hh,v 1.10 2001-07-11 09:59:31 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4BSplineCurve
//
// Class description:
// 
// Definition of a generic BSpline curve.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __BSPLINECURVE_H
#define __BSPLINECURVE_H 

#include "G4Curve.hh"
#include "G4Point3DVector.hh"

class G4ControlPoints;
class G4KnotVector;

class G4BSplineCurve : public G4Curve
{
public:

  typedef G4std::vector<G4double> G4doubleVector;

public:  // with description

  G4BSplineCurve();
    // Default constructor, must be followed by call to Init(...) method.

  virtual ~G4BSplineCurve();
    // Virtual destructor.

  G4BSplineCurve(const G4BSplineCurve& right);
  G4BSplineCurve& operator=(const G4BSplineCurve& right);
    // Copy constructor and assignment operator.

  virtual G4Curve* Project(const G4Transform3D& tr=
			   G4Transform3D::Identity);
    // Transforms and projects all control points.

  virtual G4double  GetPMax() const;
  virtual G4Point3D GetPoint(G4double param) const;
  virtual G4double  GetPPoint(const G4Point3D& p) const;
    // Accessors methods.

  void Init(G4int degree0, G4Point3DVector* controlPointsList0,
	    G4doubleVector* knots0, G4doubleVector* weightsData0);
  G4int GetDegree() const;
  const G4Point3DVector* GetControlPointsList() const;
  const G4doubleVector*  GetKnots() const;
  const G4doubleVector*  GetWeightsData() const;
    // Get/Set methods for the geometric data.
    // The "knots" vector contains each knot multiplicity Times.
    // The object is responsible for deleting the containers passed
    // to the Init method.

public:  // without description

  virtual G4bool Tangent(G4CurvePoint& cp, G4Vector3D& v);
    // Returns false. Empty implementation.

  virtual G4int IntersectRay2D(const G4Ray& ray);
    // Returns 0. Empty implementation.

  //virtual void IntersectRay2D(const G4Ray& ray, G4CurveRayIntersection& is);

protected:

  virtual void InitBounded();

  //public:
  //void ProjectCurve(const G4Plane&, const G4Plane&);   
  //int Inside(const G4Point3d&, const G4Ray&);
  //void CalcCurvePlaneNormal();

protected:

  // geometric data:
  // knots - contains each knot multiplicity Times.
  // knot_multiplicities is not needed, weightsData might be 0
  // curve_form, closed_curve, self_intersect is not used,
  // as they are unreliable sources of information.
  //
  G4int degree;
  G4Point3DVector* controlPointsList;
  G4doubleVector*  knots;
  G4doubleVector*  weightsData;

};

#include "G4BSplineCurve.icc"

#endif
