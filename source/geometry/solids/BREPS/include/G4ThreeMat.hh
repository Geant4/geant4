// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ThreeMat.hh,v 1.4 2000-02-16 12:02:52 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/*  G4ThreeMat.h,v 1.7 1993/12/30 02:15:55 rensing Exp  */
//  File: G4ThreeMat.h
//  Author: Alan Breakstone

// Contents ---------------------------------------------------------------
//
//	G4ThreeMat
//
//  Description
//
//      Defines the class G4ThreeMat for three by three matrices
//
//
// End --------------------------------------------------------------------


// Interface Dependencies -------------------------------------------------


#ifndef __THREEMAT_H
#define __THREEMAT_H

#include "G4Vector3D.hh"

// End Interface Dependencies -------------------------------------------

// Class //

class G4ThreeMat
{
//  The elements exist individually and are also aggregated into
//  rows and columns to use operations already written for the G4Vector3Dc 
//  class.
	G4double element[3][3];
	G4Vector3D row[3], column[3];
public:
//  default constructor
	G4ThreeMat();
//  Normal constructors with a 3 x 3 arG4Ray argument
	G4ThreeMat( G4double a[3][3] );
//  destructor
	virtual ~G4ThreeMat() {};
//  copy constructor
	G4ThreeMat( const G4ThreeMat& m );
//  function to return class name
	virtual const char* NameOf() const { return "G4ThreeMat"; }
//  printing functions (derived classes do not need to overwrite operator <<)
	friend G4std::ostream& operator<<( G4std::ostream& os, const G4ThreeMat& m );
	virtual void PrintOn( G4std::ostream& os = G4cout ) const;
//  equality operator
	int operator==( const G4ThreeMat& m );
//
//  overload operators =, +, -, +=, -=, *
//
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
  //
  //  Determinant of matrix
  G4double Determinant();
};

#endif
