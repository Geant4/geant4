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
// $Id$
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
	    std::vector<G4double>* knots0, std::vector<G4double>* weightsData0);
  G4int GetDegree() const;
  const G4Point3DVector* GetControlPointsList() const;
  const std::vector<G4double>*  GetKnots() const;
  const std::vector<G4double>*  GetWeightsData() const;
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
  std::vector<G4double>*  knots;
  std::vector<G4double>*  weightsData;

};

#include "G4BSplineCurve.icc"

#endif
