// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ThreeMat.hh,v 1.5 2000-08-28 08:57:50 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ThreeMat
//
// Class description:
// 
// Defines the class G4ThreeMat for three by three matrices.

// Author: Alan Breakstone
// Adaptation by: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __THREEMAT_H
#define __THREEMAT_H

#include "G4Vector3D.hh"

class G4ThreeMat
{

public:  // with description

  G4ThreeMat();
    // Default constructor.

  G4ThreeMat( G4double a[3][3] );
    // Normal constructors with a 3 x 3 array argument.

  virtual ~G4ThreeMat() {};
    // Destructor.

  G4ThreeMat( const G4ThreeMat& m );
    // Copy constructor.

  virtual const char* NameOf() const;
    // Returns the class name.

  friend G4std::ostream& operator<<( G4std::ostream& os, const G4ThreeMat& m );
  virtual void PrintOn( G4std::ostream& os = G4cout ) const;
    // Printing functions (derived classes do not need to overwrite operator <<).

  G4int operator==( const G4ThreeMat& m );
    // Equality operator.

  void operator=( const G4ThreeMat& m );
  G4ThreeMat operator+() const { return *this; };
  G4ThreeMat operator-();
  G4ThreeMat operator+=( const G4ThreeMat& m2 );
  G4ThreeMat operator-=( const G4ThreeMat& m2 );
  friend G4ThreeMat operator+( const G4ThreeMat& m1, const G4ThreeMat& m2 );
  friend G4ThreeMat operator-( const G4ThreeMat& m1, const G4ThreeMat& m2 );
  friend G4ThreeMat operator*( G4double x, const G4ThreeMat& m );
  friend G4Vector3D operator*( const G4ThreeMat& m, const G4Vector3D& v );
  friend G4ThreeMat operator*( const G4ThreeMat& m1, const G4ThreeMat& m2 );
    // Overload operators =, +, -, +=, -=, *

  G4double Determinant();
    // Determinant of matrix.

private:

  G4double element[3][3];
  G4Vector3D row[3], column[3];
    // The elements exist individually and are also aggregated into rows and
    // columns to use operations already written for the G4Vector3Dc class.
};

#endif
