// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OsloMatrix.cc,v 1.2 1999-12-15 14:50:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4OsloMatrix.hh"


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

