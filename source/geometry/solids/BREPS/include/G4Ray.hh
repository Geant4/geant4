// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Ray.hh,v 1.3 2000-08-28 08:57:48 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4Ray
//
// Class description:
// 
// Definition of a generic ray.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4Ray_h
#define __G4Ray_h 1

#include "G4Point3D.hh"
#include "G4PointRat.hh"
#include "G4Vector3D.hh"
#include "G4Plane.hh"


class G4Ray 
{
  
public:  // with description

  G4Ray();
  G4Ray(const G4Point3D& start0, const G4Vector3D& dir0);
  ~G4Ray();
    // Contructors and destructor.

  inline G4Point3D GetPoint(G4double i) const;
  inline G4double GetPPoint(const G4Point3D& p) const;
  inline const G4Vector3D& GetDir() const;
  inline const G4Point3D& GetStart() const;
  inline void SetDir(const G4Vector3D& dir0);
  inline void SetStart(const G4Point3D& start0);
  const G4Plane& GetPlane(G4int number_of_plane) const;
    // Get/Set methods of geometrical data.

  void RayCheck();
    // Makes sure that the vector has unit length.

  void CreatePlanes();
    // Creates two orthogonal planes (plane1,plane2), the ray (rray) 
    // being situated in the intersection of the planes. The planes are 
    // used to project the surface (nurb) in two dimensions.

  static G4int CalcPlane3Pts( G4Plane& plane, const G4Point3D& a,
			      const G4Point3D& b, const G4Point3D& c );
    // Finds the equation of a G4Plane that contains three points.
    // Note that Normal vector created is expected to point outside.
    // This follows the outward-pointing Normal convention, and the
    // right-hand rule for cross products.
    //
    //
    //                      C
    //                      *
    //                      |\
    //                      | \
    //         ^     N      |  \
    //         |      \     |   \
    //         |       \    |    \
    //         |C-A     \   |     \
    //         |         \  |      \
    //         |          \ |       \
    //                     \|        \
    //                      *---------*
    //                      A         B
    //                         ----->
    //                          B-A
    //
    // If the points are given in the order A B C (eg, *counter*-clockwise),
    // then the outward pointing surface Normal N = (B-A) x (C-A).
    //
    // Explicit return value:
    //       0      OK
    //      -1      Failure.  At least two of the points were not distinct,
    //              or all three were colinear.
    //
    // Implicit return argument:
    //      plane   The G4Plane equation coefficients are stored here.

  inline G4double P2(G4double x);
  inline G4int NearZero(G4double val, G4double epsilon);
  void MatVecOrtho(register G4Vector3D &out, register const G4Vector3D in);
    // Utility methods.

  inline void Vsetall(G4Vector3D &a, G4double s);
    // Sets all elements of vector to the same scalar value.
  
  static inline void Vcross(G4Plane &a, const G4Vector3D &b,
                                        const G4Vector3D &c);
    // Cross product of 'b' and 'c'. Stores result in 'a' (G4Plane).

  static inline void Vcross(G4Vector3D &a, const G4Vector3D &b,
                                           const G4Vector3D &c);
    // Cross product of 'b' and 'c'. Stores result in 'a' (G4Vector3D).
  
  static inline void Vmove(G4Point3D &a, const G4Point3D &b);
    // Sets 'a' equal to 'b'.

  static inline void Vadd2(G4Point3D &a, const G4Point3D &b,
                                         const G4Vector3D &c );
    // Adds vector 'c' to 'b'. Stores result in 'a'.
  
  static inline void Vsub2(G4Vector3D &a, const G4Point3D &b,
                                          const G4Point3D &c);
    // Subtracts vector 'c' from 'b'. Stores result in 'a'.

  static inline void Vscale(G4Plane& a, const G4Plane& b, G4double c);
    // Scales vector at `b' by scalar `c'. Stores result in `a'.

  static inline  G4double Vdot(const G4Plane &a, const G4Point3D &b);
    // Computes dot product of vectors at `a' and `b'.
  
  static inline G4double Magsq(const G4Plane &a);
    // Returns scalar Magnitude squared of vector at `a'.
  
  static inline G4double Magnitude(const G4Plane &a);
    // Returns scalar Magnitude of vector at `a'.

public:  // without description

  void Init(const G4Point3D& start0, const G4Vector3D& dir0);
    // Initialisation of a G4Ray (called by constructor).

private:

  G4Point3D  start;
  G4Vector3D dir;

  G4double   r_min;		// entry Dist to bounding sphere 
  G4double   r_max;		// exit Dist from bounding sphere 
  
  G4Plane    plane1, plane2; 

};

#include "G4Ray.icc"

#endif
