// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OsloMatrix.cc,v 1.3 2000-08-28 08:57:58 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4OsloMatrix.cc
//
// ----------------------------------------------------------------------

#include "G4OsloMatrix.hh"

G4OsloMatrix::G4OsloMatrix()
{
  o_vec = (G4KnotVector*)0;
  next  = (G4OsloMatrix*)0;
}

G4OsloMatrix::G4OsloMatrix(G4int vec_size, G4int offsetparam, G4int osizeparam)
{
  next   = (G4OsloMatrix*)0;
  o_vec  = new G4KnotVector(vec_size);
  offset = offsetparam;
  osize  = osizeparam;
}

G4OsloMatrix::~G4OsloMatrix()
{
  delete o_vec;
}

Matrix::Matrix()
{
  nr=nc=0;
  data=0;
}


Matrix::Matrix(int rows, int columns)
{
  nr=rows; 
  nc=columns; 
  data = new G4double[nr*nc];
  
  for(int a =0; a<nr*nc;a++) 
    data[a]=0;
}


Matrix::Matrix(G4double vec[])
{
  nr = 4;
  nc = 4; 
  data = new G4double[nr*nc];
  
  for(int a=0;a<nr*nc;a++)
    data[a]=vec[a];
}

Matrix::~Matrix(){;}

