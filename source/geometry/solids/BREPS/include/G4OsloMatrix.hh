//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4OsloMatrix.hh,v 1.5 2001-07-11 09:59:35 gunter Exp $
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

  inline G4int GetOffset() const;
  inline G4int GetSize() const;
  inline G4OsloMatrix* GetNextNode();
  inline G4KnotVector* GetKnotVector();
  inline void SetOffset(G4int);
  inline void SetSize(G4int);
  inline void SetNextNode(G4OsloMatrix*);
  inline void SetKnotVector(G4KnotVector*);
    // Accessors

private:

  G4OsloMatrix(const G4OsloMatrix&);
  G4OsloMatrix& operator=(const G4OsloMatrix&);
    // Private copy constructor and assignment operator.

private:

  G4OsloMatrix* next;
  G4int offset;
  G4int osize;
  G4KnotVector* o_vec;
};


class G4PointMatrix
{
public:

  G4PointMatrix();
  G4PointMatrix(G4int, G4int);
  G4PointMatrix(G4double[]);
    // Constructors

  ~G4PointMatrix();
    // Destructor

  inline G4int GetRows() const;
  inline G4int GetCols() const;

  inline void put(G4int i, G4int j, G4double x);
    // Puts control point into matrix location (i,j)
  
  inline G4double get(G4int i, G4int j) const;
    // Retrieves control point from matrix location (i,j)

private:

  G4PointMatrix(const G4PointMatrix&);
  G4PointMatrix& operator=(const G4PointMatrix&);
    // Private copy constructor and assignment operator.

private:				      

  G4double* data;
  G4int nr, nc;
};

#include "G4OsloMatrix.icc"

#endif
