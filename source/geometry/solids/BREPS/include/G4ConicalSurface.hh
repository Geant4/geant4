// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ConicalSurface.hh,v 1.5 2000-08-28 08:57:44 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ConicalSurface
//
// Class Description:
//   
// A G4ConicalSurface is a semi-infinite conical surface defined by 
// an axis and an opening angle, defined as the angle between the axis 
// and the conical surface, with the origin being the apex of the cone.

// The code for G4ConicalSurface has been derived from the original
// implementation in the "Gismo" package.
//
// Author: A.Breakstone
// Adaptation: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4CONICALSURFACE_H
#define __G4CONICALSURFACE_H

#include "G4Surface.hh"
 

class G4ConicalSurface : public G4Surface
{

public:  // with description

  G4ConicalSurface();
    // Default constructor:
    //   default axis is ( 1.0, 0.0, 0.0 ),
    //   default angle is 1.0 radians.

  G4ConicalSurface( const G4Point3D& o, const G4Vector3D& a, G4double e ); 
    // Normal constructor:
    //   first argument is the origin of the G4ConicalSurface
    //   second argument is the axis of the G4ConicalSurface
    //   third argument is the angle of the G4ConicalSurface.

  virtual ~G4ConicalSurface();
    // Virtual destructor.

  inline G4int operator==( const G4ConicalSurface& c );
    // Equality operator.

  inline G4String GetEntityType();
    // Returns type identifier of the shape.

  virtual const char* NameOf() const;
    // Returns the class name.

  virtual void PrintOn( G4std::ostream& os = G4cout ) const;
    // Printing function, streaming surface's attributes.

  virtual G4double HowNear( const G4Vector3D& x ) const;
    // Returns the distance from a point to a semi-infinite G4ConicalSurface.
    // The point x is the (input) argument.
    // The distance is positive if the point is Inside, negative if it
    // is outside

  void CalcBBox();
    // Computes the bounding-box.

  G4int Intersect( const G4Ray& ry );
    // Returns the distance along a Ray (straight line with G4Vector3D) to
    // leave or enter a G4ConicalSurface.
    // If the G4Vector3D of the Ray is opposite to that of the Normal to
    // the G4ConicalSurface at the intersection point, it will not leave the 
    // G4ConicalSurface.
    // Similarly, if the G4Vector3D of the Ray is along that of the Normal 
    // to the G4ConicalSurface at the intersection point, it will not enter the
    // G4ConicalSurface.
    // This method is called by all finite shapes sub-classed to 
    // G4ConicalSurface.
    // A negative result means no intersection.
    // If no valid intersection point is found, set the distance
    // and intersection point to large numbers.

  virtual G4Vector3D SurfaceNormal( const G4Point3D& p ) const;	
    // Returns the Normal unit vector to the G4ConicalSurface at a point p 
    // on (or nearly on) the G4ConicalSurface.

  virtual G4int Inside( const G4Vector3D& x ) const;
    // Returns 1 if the point x is Inside the G4ConicalSurface, 0 otherwise.
    // Outside means that the distance to the G4ConicalSurface would be
    // negative. Uses the HowNear() function to calculate this distance.

  virtual G4int WithinBoundary( const G4Vector3D& x ) const;
    // Returns 1 if point x is on the G4ConicalSurface, otherwise return zero
    // Since a G4ConicalSurface is infinite in extent, the function 
    // will just check if the point is on the G4ConicalSurface (to the surface
    // precision).

  virtual G4double Scale() const;
    // Function overwritten by finite-sized derived classes which returns
    // a radius, unless it is zero, in which case it returns the smallest
    // non-zero dimension. 
    // Since a semi-infinite cone has no Scale associated with it, it returns
    // the arbitrary number 1.0.
    // Used for Scale-invariant tests of surface thickness.

  inline G4Vector3D GetAxis() const;
  inline G4double GetAngle() const;
    // Return the axis and angle of the G4ConicalSurface.

  void SetAngle( G4double e );
    // Changes the angle of the G4ConicalSurface.
    // Requires angle to range from 0 to PI/2.

public:  // without description

/*
  G4ConicalSurface( const G4ConicalSurface& c );
    // Copy constructor.

  virtual G4double distanceAlongRay( G4int which_way, const G4Ray* ry,
                                     G4Vector3D& p ) const;
    // Returns the distance along a Ray to enter or leave a G4ConicalSurface.  
    //	 The first (input) argument is +1 to leave or -1 to enter
    //   The second (input) argument is a pointer to the Ray
    //   The third (output) argument returns the intersection point.

  virtual G4double distanceAlongHelix( G4int which_way, const Helix* hx,
                                       G4Vector3D& p ) const;
    // Returns the distance along a Helix to enter or leave a G4ConicalSurface.  
    //   The first (input) argument is +1 to leave or -1 to enter
    //   The second (input) argument is a pointer to the Helix
    //   The third (output) argument returns the intersection point.

  G4Vector3D Normal( const G4Vector3D& p ) const;
    // Returns the Normal unit vector to a G4ConicalSurface 
    // at a point p on (or nearly on) the G4ConicalSurface.

  virtual void rotate( G4double alpha, G4double beta, 
                       G4double gamma, G4ThreeMat& m, G4int inverse ); 
    // Rotates the G4ConicalSurface (angles are assumed to be given in
    // radians), arguments:
    // - first about global x-axis by angle alpha,
    // - second about global y-axis by angle beta,
    // - third about global z-axis by angle gamma,
    // - fourth (output) argument gives the calculated rotation matrix,
    // - fifth (input) argument is an integer flag which if
    //   non-zero reverses the order of the rotations.

  virtual void rotate( G4double alpha, G4double beta, 
                       G4double gamma, G4int inverse ); 
    // Rotates the G4ConicalSurface (angles are assumed to be given in
    // radians), arguments: 
    // - first about global x-axis by angle alpha,
    // - second about global y-axis by angle beta,
    // - third about global z-axis by angle gamma,
    // - fourth (input) argument is an integer flag which if
    //   non-zero reverses the order of the rotations.

private:

  virtual G4double gropeAlongHelix( const Helix* hx ) const;
    // Private function to use a crude technique to find the intersection
    // of a Helix with a G4ConicalSurface.  It returns the turning angle
    // along the Helix at which the intersection occurs or -1.0 if no
    // intersection point is found.
    // The argument to the call is the pointer to the Helix.
*/

private:

  G4Vector3D axis;
    // Direction of axis of G4ConicalSurface (unit vector).

  G4double   angle;
    // Half opening angle of G4ConicalSurface, in radians
    // range is 0 < angle < PI/2.

};

#include "G4ConicalSurface.icc"

#endif
