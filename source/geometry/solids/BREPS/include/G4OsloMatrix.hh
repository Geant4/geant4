// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OsloMatrix.hh,v 1.3 2000-08-28 08:57:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4OsloMatrix
//
// Class description:
// 
// Utility class for the definition of a matrix of knots.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4OsloMatrix_h
#define __G4OsloMatrix_h 1

#include "G4KnotVector.hh"


class G4OsloMatrix 
{

public:  // with description

  G4OsloMatrix();
  G4OsloMatrix(G4int vec_size, G4int offsetparam, G4int osizeparam);
  ~G4OsloMatrix();
    // Constructors & destructor

public:  // without description

  G4OsloMatrix * next;
  G4int offset;
  G4int osize;
  G4KnotVector *o_vec;
};


class Matrix
{
public:

  Matrix();
  Matrix(G4int, G4int);
  Matrix(G4double[]);
    // Constructors

  ~Matrix();
    // Destructor

  inline G4int GetRows() const;
  inline G4int GetCols() const;

  inline void put(G4int i, G4int j, G4double x);
    // Puts control point into matrix location (i,j)
  
  inline G4double get(G4int i, G4int j) const;
    // Retrieves control point from matrix location (i,j)


private:				      

  G4double* data;
  G4int nr, nc;
};

#include "G4OsloMatrix.icc"

#endif
