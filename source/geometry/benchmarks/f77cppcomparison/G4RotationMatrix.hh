// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RotationMatrix.hh,v 1.2 1999-12-15 14:49:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Rotation Matrix class, converted from CLHEP:
// Author: Leif Lonnblad

// History:
// 30.11.94 P.Kent: Added phiX/thetaX (etc) + IsIdentity functions

#ifndef G4ROTATIONMATRIX_HH
#define G4ROTATIONMATRIX_HH

#include "globals.hh"

//#include "G4ThreeVector.hh"
class G4ThreeVector;

class G4RotationMatrix {

public:
    
    // Default constructor. Gives a unit matrix.
    inline G4RotationMatrix()
	: xx(1.0), xy(0.0), xz(0.0),
	  yx(0.0), yy(1.0), yz(0.0),
	  zx(0.0), zy(0.0), zz(1.0)
    {;}
    
    // Copy constructor.
    inline G4RotationMatrix(const G4RotationMatrix & m)
	: xx(m.xx), xy(m.xy), xz(m.xz),
	  yx(m.yx), yy(m.yy), yz(m.yz),
	  zx(m.zx), zy(m.zy), zz(m.zz)
    {;}


    inline G4RotationMatrix & operator = (const G4RotationMatrix & m);
    // Assignment.

    inline G4bool operator == (const G4RotationMatrix &m) const;
    // Comparison

    inline G4ThreeVector operator * (const G4ThreeVector &) const;
    // Multiplication with a d3Vector

    G4RotationMatrix operator * (const G4RotationMatrix &) const;
    inline G4RotationMatrix & operator *= (const G4RotationMatrix &);
    inline G4RotationMatrix & transform(const G4RotationMatrix &);
    // Matrix multiplication.
    // Note a *= b; <=> a = a * b; while a.transform(b); <=> a = b * a;

    inline G4RotationMatrix inverse() const;
    // Returns the inverse.

    inline G4RotationMatrix & invert();
    // Inverts the Rotation matrix
    
    G4RotationMatrix & rotateX(double);
    // Rotation around the x-axis.

    G4RotationMatrix & rotateY(double);
    // Rotation around the y-axis.
    
    G4RotationMatrix & rotateZ(double);
    // Rotation around the z-axis.
    
    G4RotationMatrix & rotate(double angle, const G4ThreeVector & axis);
    inline G4RotationMatrix & rotate(double angle, const G4ThreeVector * axis);
    // Rotation around a specified vector.

// Function to return angles (RADS) made by rotated axes against original axes
    inline double phiX() const;
    inline double phiY() const;
    inline double phiZ() const;

    inline double thetaX() const;
    inline double thetaY() const;
    inline double thetaZ() const;

    inline G4bool isIdentity() const;
// IsIdentity function returns true if identity  matrix


protected:

  inline G4RotationMatrix(double, double, double, double, double,
		     double, double, double, double);
  // Protected constructor;.

  G4double xx, xy, xz, yx, yy, yz, zx, zy, zz;
  // The matrix elements.

};

#include "G4ThreeVector.hh"

#include "G4RotationMatrix.icc"

#endif





