// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OsloMatrix.cc,v 1.4 2000-11-08 14:22:10 gcosmo Exp $
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

G4PointMatrix::G4PointMatrix()
{
  nr=nc=0;
  data=0;
}


G4PointMatrix::G4PointMatrix(int rows, int columns)
{
  nr=rows; 
  nc=columns; 
  data = new G4double[nr*nc];
  
  for(int a =0; a<nr*nc;a++) 
    data[a]=0;
}


G4PointMatrix::G4PointMatrix(G4double vec[])
{
  nr = 4;
  nc = 4; 
  data = new G4double[nr*nc];
  
  for(int a=0;a<nr*nc;a++)
    data[a]=vec[a];
}

G4PointMatrix::~G4PointMatrix(){;}

