// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Line.hh,v 1.1 1999-01-07 16:07:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef __LINE_H
#define __LINE_H 

#include "G4Curve.hh"

class G4Line : public G4Curve
{
public:

  G4Line ();
  virtual ~G4Line ();

  virtual G4Curve* Project(const G4Transform3D& tr = G4Transform3D::Identity);

  virtual G4bool Tangent(G4CurvePoint& cp, G4Vector3D& vec);

  virtual void IntersectRay2D(const G4Ray& ray, G4CurveRayIntersection& is);

  virtual G4double  GetPMax();
  virtual G4Point3D GetPoint(G4double param);
  virtual G4double  GetPPoint(const G4Point3D& pt);

  // Get/Set for the geometric data
  void Init(const G4Point3D& pnt0, const G4Vector3D& dir0);
  G4Point3D  GetPnt() const;
  G4Vector3D GetDir() const;


protected:

  virtual void InitBounded();

 
private:
  // For the Inside function
  //inline int Sign(G4double a, G4double b){return((a>=0&&b>=0)||(a<0&&b<0));}
  
  // geometric data
  G4Point3D  pnt;
  G4Vector3D dir;
  G4Vector3D invDir;   // dir / |dir|^2 always
  G4Vector3D v;        // dir / |dir| always

};

#include "G4Line.icc"

#endif
