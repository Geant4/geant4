// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Ray.hh,v 1.2 1999-12-15 14:49:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef __G4Ray_h
#define __G4Ray_h 1

#include "G4Point3D.hh"
#include "G4PointRat.hh"
#include "G4Vector3D.hh"
#include "G4Plane.hh"


class G4Ray 
{
  
public:

  G4Ray();

  G4Ray(const G4Point3D& start0, const G4Vector3D& dir0);

  void Init(const G4Point3D& start0, const G4Vector3D& dir0);

  G4Point3D GetPoint(G4double i) const;

  G4double GetPPoint(const G4Point3D& p) const;

  const G4Vector3D& GetDir() const;

  const G4Point3D& GetStart() const;

  void SetDir(const G4Vector3D& dir0);

  void SetStart(const G4Point3D& start0);

private:

  G4Point3D  start;
  G4Vector3D dir;

  G4double   r_min;		// entry Dist to bounding sphere 
  G4double   r_max;		// exit Dist from bounding sphere 
  
  G4Plane    plane1, plane2; 


public:

  const G4Plane& GetPlane(const int number_of_plane)const;//1 or 2

  void RayCheck();

  void CreatePlanes();

  static int CalcPlane3Pts( G4Plane &plane1, const G4Point3D& a,
			    const G4Point3D& b, const G4Point3D& c );

  inline G4double P2(const G4double x) {return(x*x);};

  void MatVecOrtho( register G4Vector3D &out,  register const G4Vector3D in );

  inline int NearZero(const G4double val, const G4double epsilon) 
  {
    return ( ((val) > -epsilon) && ((val) < epsilon) );
  }

  static inline void Vcross(G4Plane &a, 
			    const G4Vector3D &b, const G4Vector3D &c) 
  { 
    a.a = b.y()  * c.z()  - b.z()  * c.y() ;
    a.b = b.z()  * c.x()  - b.x()  * c.z() ;
    a.c = b.x()  * c.y()  - b.y()  * c.x() ;
  }

  static inline void Vcross(G4Vector3D &a, 
			    const G4Vector3D &b, const G4Vector3D &c) 
  { 
    a.setX(b.y()  * c.z()  - b.z()  * c.y()) ;
    a.setY(b.z()  * c.x()  - b.x()  * c.z()) ;
    a.setZ(b.x()  * c.y()  - b.y()  * c.x()) ;
  }
  
  inline void Vmove(G4Point3D &a, const G4Point3D &b) 
  { 
    a.setX(b.x());
    a.setY(b.y());
    a.setZ(b.z());
  }

  inline void Vadd2(G4Point3D &a, const G4Point3D &b, const G4Vector3D &c ) 
  {
    a.setX(b.x() + c.x()) ;
    a.setY(b.y() + c.y()) ;
    a.setZ(b.z() + c.z()) ;
  }	
  
  static inline void Vsub2(G4Vector3D &a, 
			   const G4Point3D &b, const G4Point3D &c) 
  {
    a.setX(b.x() - c.x());
    a.setY(b.y() - c.y());
    a.setZ(b.z() - c.z());
  }

  // Set all elements of vector to same scalar value
  inline void Vsetall(G4Vector3D &a, G4double s) 
  {
    a.setX(s); a.setY(s); a.setZ(s);
  }
  
  // Scale vector at `b' by scalar `c', Store result at `a'
  static inline void Vscale(G4Plane& a, const G4Plane& b, const G4double c) 
  { 
    a.a = b.a * c;
    a.b = b.b * c;
    a.c = b.c * c;
  }

  // Compute dot product of vectors at `a' and `b' 
  static inline  G4double Vdot(const G4Plane &a, const G4Point3D &b) 
  {
    return (a.a * b.x() + 
	    a.b * b.y() + 
	    a.c * b.z());
  }
  
  // Return scalar Magnitude squared of vector at `a'
  static inline G4double Magsq(const G4Plane &a) 
  {
    return ( a.a * a.a + a.b * a.b + a.c *a.c );
  }
  
  // Return scalar Magnitude of vector at `a' 
  static inline G4double Magnitude(const G4Plane &a) 
  {
    return (sqrt( Magsq( a )) );
  }

};

#include "G4Ray.icc"

#endif
















