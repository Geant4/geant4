// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OsloMatrix.cc,v 1.1 1999-01-07 16:07:44 gunter Exp $
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

