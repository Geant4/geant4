// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SphericalSurface.hh,v 1.7 2000-11-08 14:22:04 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4SphericalSurface
//
// Class description:
// 
// Definition of a spherical surface.

// The code for G4SphericalSurface has been derived from the original
// implementation in the "Gismo" package.
//
// Authors: L.Lim, A.Breakstone.
// Adaptation: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4SpheShell_H
#define __G4SpheShell_H

#include "G4Surface.hh"
#include "G4ThreeMat.hh"

class G4SphericalSurface : public G4Surface
{

public:  // with description

  G4SphericalSurface();
    // Default constructor.

  G4SphericalSurface( const G4Vector3D& o, 
		      const G4Vector3D& xhat, const G4Vector3D& zhat,
		      G4double r, 
		      G4double ph1, G4double ph2,
		      G4double th1, G4double th2 ); 
    // Normal constructor:
    //   first argument is the origin of the G4SphericalSurface
    //   second argument is the axis of the G4SphericalSurface
    //		which defines azimuthal angle equals zero
    //   third argument is the axis of the G4SphericalSurface
    //          which defines polar angle equals zero
    //   fourth argument is the radius of the G4SphericalSurface
    //   fifth argument is the lower azimuthal angle limit of the surface
    //   sixth argument is the upper azimuthal angle limit of the surface
    //   seventh argument is the lower polar angle limit of the surface
    //   eigth argument is the upper polar angle limit of the surface

  virtual ~G4SphericalSurface();
    // Destructor.

  inline G4int operator==( const G4SphericalSurface& s );
    // Equality operator.

  inline G4String GetEntityType() const;
    // Returns the type identifier.

  virtual const char* NameOf() const;
    // Returns the class name.

  virtual void PrintOn( G4std::ostream& os = G4cout ) const;
    // Printing function, streaming surface's attributes.

  G4int Intersect(const G4Ray&);
    // Returns the distance along a Ray (straight line with G4Vector3D) to
    // leave or enter a G4SphericalSurface.
    // If the G4Vector3D of the Ray is opposite to that of the Normal to
    // the G4SphericalSurface at the intersection point, it will not leave the
    // G4SphericalSurface.
    // Similarly, if the G4Vector3D of the Ray is along that of the Normal 
    // to the G4SphericalSurface at the intersection point, it will not enter 
    // the G4SphericalSurface.
    // This method is called by all finite shapes sub-classed to 
    // G4SphericalSurface.
    // A negative result means no intersection.
    // If no valid intersection point is found, set the distance
    // and intersection point to large numbers.

  void CalcBBox();
    // Computes the bounding-box.

  inline void Comp(G4Vector3D& v, G4Point3D& min , G4Point3D& max);
    // Compares the x,y and z values of v and min
    // versus v and max. min/max-values are replaced if 
    // greater/smaller than v-values.

  virtual G4double HowNear( const G4Vector3D& x ) const;
    // Returns the distance from a point to a G4SphericalSurface
    // The point x is the (input) argument.
    // The distance is positive if the point is Inside, negative if it
    // is outside
  
  virtual G4Vector3D SurfaceNormal( const G4Point3D& p ) const;
    // Returns the Normal unit vector to the G4SphericalSurface at a point p 
    // on (or nearly on) the G4SphericalSurface.

  virtual G4int Inside( const G4Vector3D& x ) const;
    // Returns 1 if the point x is Inside the G4SphericalSurface, 0 otherwise.

  virtual G4int WithinBoundary( const G4Vector3D& x ) const;
    // Returns 1 if the point x is within the boundary, 0 otherwise.

  virtual G4double Scale() const;
    // Returns the radius, unless it is zero, in which case it
    // returns 1.  Used for Scale-invariant tests of surface thickness.

  virtual G4double Area() const; 
    // Calculates the area of a G4SphericalSurface.

  virtual void resize( G4double r, G4double ph1, G4double ph2, 
		       G4double th1, G4double th2);
    // Resizes the G4SphericalSurface to new radius and angle limits.
    //   first argument is the radius of the G4SphericalSurface
    //   second argument is the lower azimuthal angle limit of the surface
    //   third argument is the upper azimuthal angle limit of the surface
    //   fourth argument is the lower polar angle limit of the surface
    //   fifth argument is the upper polar angle limit of the surface

  inline G4Vector3D GetXAxis() const;
  inline G4Vector3D GetZAxis() const;
  inline G4double   GetRadius() const;
  inline G4double   GetPhi1() const;
  inline G4double   GetPhi2() const;
  inline G4double   GetTheta1() const;
  inline G4double   GetTheta2() const;
    // Accessors methodss to return the axes, radius, and angles
    // of the G4SphericalSurface.

public:  // without description

  virtual G4Vector3D Normal( const G4Vector3D& p ) const;
    // Returns the Normal unit vector as for SurfaceNormal().

/*
  virtual G4double distanceAlongRay( G4int which_way, const G4Ray* ry,
                                     G4ThreeVec& p ) const;
    // Returns the distance along a Ray to enter or leave a G4SphericalSurface.  
    //   The first (input) argument is +1 to leave or -1 to enter
    //   The second (input) argument is a pointer to the Ray
    //   The third (output) argument returns the intersection point.

  virtual G4double distanceAlongHelix( G4int which_way, const Helix* hx,
                                       G4ThreeVec& p ) const;
    // Returns the distance along a Helix to enter or leave a G4SphericalSurface.  
    //   The first (input) argument is +1 to leave or -1 to enter
    //   The second (input) argument is a pointer to the Helix
    //   The third (output) argument returns the intersection point.

  virtual G4Vector3D Normal( const G4Point3D& p ) const;
    // Returns the Normal unit vector to a G4SphericalSurface at a point p
    // on (or nearly on) the G4SphericalSurface.

  virtual void rotate( G4double alpha, G4double beta, 
                       G4double gamma, G4ThreeMat& m, G4int inverse );
    // Rotates the G4SphericalSurface (angles are assumed to be given in
    // radians), arguments:
    // - first about global x_axis by angle alpha,
    // - second about global y-axis by angle beta,
    //	- third about global z_axis by angle gamma,
    //	- fourth (output) argument gives the calculated rotation matrix,
    // - fifth (input) argument is an integer flag which if
    //   non-zero reverses the order of the rotations.

  virtual void rotate( G4double alpha, G4double beta, 
                       G4double gamma, G4int inverse ); 
    // Rotates the G4SphericalSurface (angles are assumed to be given in
    // radians), arguments:
    //	- first about global x_axis by angle alpha,
    // - second about global y-axis by angle beta,
    // - third about global z_axis by angle gamma,
    // - fourth (input) argument is an integer flag which if
    //   non-zero reverses the order of the rotations.
*/

protected:  // with description

  G4Vector3D x_axis;
    // Direction (unit vector) of axis of G4SphericalSurface
    // which defines azimuthal angle of zero.

  G4Vector3D z_axis;
    // Direction (unit vector) of axis of G4SphericalSurface
    // which defines polar angle of zero.
	
  G4double radius;
    // Radius of G4SphericalSurface.

  G4double phi_1;
    // Lower azimuthal angle limit of G4SphericalSurface
    // (in radians).  Allowed range: 0 <= phi_1 < 2*PI.
  
  G4double phi_2;
    // Upper azimuthal angle limit of G4SphericalSurface
    // (in radians).  Allowed range: phi_1 < phi_2 <= phi_1 + 2*PI

  G4double theta_1;
    // Lower polar angle limit of G4SphericalSurface
    // (in radians).  Allowed range: 0 <= theta_1 < PI.
	
  G4double theta_2;
    // Upper polar angle limit of G4SphericalSurface
    // (in radians).  Allowed range: theta_1 < theta_2 <= theta_1 + PI.

private:

  G4SphericalSurface(const G4SphericalSurface&);
  G4SphericalSurface& operator=(const G4SphericalSurface&);
    // Private copy constructor and assignment operator.

  // virtual G4double gropeAlongHelix( const Helix* hx ) const;
    // Private function to use a crude technique to find the intersection
    // of a Helix with a G4SphericalSurface.  It returns the turning angle 
    // along the Helix at which the intersection occurs or -1.0 if no
    // intersection point is found.  The argument to the call is the pointer
    // to the Helix.

};

#include "G4SphericalSurface.icc"

#endif
