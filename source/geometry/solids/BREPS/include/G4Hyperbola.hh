// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Hyperbola.hh,v 1.1 1999-01-07 16:07:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef __HYPERBOLICCURVE_H
#define __HYPERBOLICCURVE_H 

#include "G4Conic.hh"

class G4Hyperbola : public G4Conic
{
public:

  G4Hyperbola();
  ~G4Hyperbola();
  
  virtual G4Curve* Project(const G4Transform3D& tr= G4Transform3D::Identity);

  virtual G4bool Tangent(G4CurvePoint& cp, G4Vector3D& v);

  virtual void IntersectRay2D(const G4Ray& ray, G4CurveRayIntersection& is);

  virtual G4double  GetPMax();
  virtual G4Point3D GetPoint(G4double param);
  virtual G4double  GetPPoint(const G4Point3D& p);

  // STEP
  G4Hyperbola(STEPentity& Ent);
  G4Hyperbola(STEPentity& Ent, InstMgr&);        

  //G4Hyperbola(G4Point3d, G4Point3d, G4Point3d, 
  //            G4Point3d,G4double, G4double );
  //G4Point3d EvaluateByParameterValue(const G4double u);
  //G4Point3d GetBoundMax();
  //G4Point3d GetBoundMin();    

  // Get/Set for the geometric data
  void Init(G4Axis2Placement3D position0,
	    G4double semiAxis0, G4double semiImagAxis0);
  
  G4double GetSemiAxis() const;
  G4double GetSemiImagAxis() const;


protected:

  virtual void InitBounded();


private:

  int Inside(const G4Point3D&, const G4Ray&);
/* L. Broglia
  G4Point3d Focus1;
  G4Point3d Focus2;    
  G4Point2d ProjFocus1;
  G4Point2d ProjFocus2;  
*/

  G4Point3D Focus1;
  G4Point3D Focus2;    
  G4Point3D ProjFocus1;
  G4Point3D ProjFocus2; 

  // geometric data
  G4double semiAxis;
  G4double semiImagAxis;

  G4double ratioAxisImagAxis;

  G4Transform3D toUnitHyperbola;

  G4double forTangent; // R_1^2/R_2^2
};

#include "G4Hyperbola.icc"

#endif
