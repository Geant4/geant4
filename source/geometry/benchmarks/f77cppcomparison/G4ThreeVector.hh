// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ThreeVector.hh,v 1.1 1999-01-08 16:31:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// ThreeVector class, converted from CLHEP:
// Authors: Leif Lonnblad and Anders Nilsson.
//
// History:
// 30.06.95 P.Kent

#ifndef G4THREEVECTOR_HH
#define G4THREEVECTOR_HH

#include "globals.hh"
#include "geomdefs.hh"

class ostream;

class G4ThreeVector {

friend class G4RotationMatrix;

public:

    // The Constructors
    inline G4ThreeVector(G4double x = 0.0, G4double y = 0.0, G4double z = 0.0)
	: dx(x), dy(y), dz(z)
    {;}
    
    inline G4ThreeVector(const G4ThreeVector & p)
	: dx(p.dx), dy(p.dy), dz(p.dz)
    {;}
    
    // the x, y and z components.
    inline G4double x() const
    {
	return dx;
    }

    inline G4double y() const
    {
	return dy;
    }

    inline G4double z() const
    {
	return dz;
    }

// Assignment.
    inline G4ThreeVector & operator = (const G4ThreeVector & p)
    {
	dx = p.x();
	dy = p.y();
	dz = p.z();
	return *this;
    }

  inline G4bool operator==(const G4ThreeVector& v) const;
    // Test for equality

  inline G4ThreeVector & operator += (const G4ThreeVector &);
  // Addition.

  inline G4ThreeVector & operator -= (const G4ThreeVector &);
  // Subtraction operator.

  inline G4ThreeVector operator - () const;
  // Unary minus.

  inline G4ThreeVector & operator *= (G4double);
  // Operators for scaling with real numbers.

  inline G4double dot(const G4ThreeVector &) const;
  // Scalar product.

  inline G4ThreeVector cross(const G4ThreeVector &) const;
  // Cross product.

  inline G4double operator () (const EAxis) const;
// Component access 


  inline G4ThreeVector unit() const;
  // the unit vector parallel to this

  inline G4double mag2() const;
  // the Magnitude squared.

  inline G4double mag() const;
  // the magnutude.

  inline G4double perp2() const;
  // The transverse component squared.

  inline G4double perp() const;
  // The transverse component.

  inline G4double perp2(const G4ThreeVector &) const;
  // The transverse component wrt. given axis squared.

  inline G4double perp(const G4ThreeVector &) const;
  // The transverse component wrt. given axis.

  inline G4double phi() const;
  // The azimuth angle.

  inline G4double theta() const;
  // The polar angle.

  inline G4double cosTheta() const;
  // Cosine of the polar angle.

  inline G4double angle(const G4ThreeVector &) const;
  // The angle w.r.t. another 3-vector.

  inline G4ThreeVector & operator *= (const G4RotationMatrix &);
  inline G4ThreeVector & transform(const G4RotationMatrix &);
  // Transformation with a Rotation matrix.

  void rotateX(G4double);
  // Rotates the G4ThreeVector around the x-axis.

  void rotateY(G4double);
  // Rotates the G4ThreeVector around the y-axis.

  void rotateZ(G4double);
  // Rotates the G4ThreeVector around the z-axis.

  void rotate(G4double, const G4ThreeVector &);
  // Rotates around the axis specified by another G4ThreeVector.


protected:

  inline void setX(G4double);
  inline void setY(G4double);
  inline void setZ(G4double);
  // Set the x, y and z components.

private:

  G4double dx, dy, dz;
  // The components.

};

ostream & operator << (ostream &, const G4ThreeVector &);
// output to a stream


#include "G4ThreeVector.icc"
// Inline functions

#endif
