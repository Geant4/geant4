// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OsloMatrix.hh,v 1.2 1999-12-15 14:49:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef __G4OsloMatrix_h
#define __G4OsloMatrix_h 1

#include "G4KnotVector.hh"


class G4OsloMatrix 
{
public:

  G4OsloMatrix()
  {
    o_vec = (G4KnotVector*)0;
    next  = (G4OsloMatrix*)0;
  };
  
  G4OsloMatrix(int vec_size, int offsetparam, int osizeparam)
  {
    next   = (G4OsloMatrix*)0;
    o_vec  = new G4KnotVector(vec_size);
    offset = offsetparam;
    osize  = osizeparam;
  }
	
  ~G4OsloMatrix() { delete o_vec; }
	
  G4OsloMatrix * next;
  int offset;
  int osize;
  G4KnotVector *o_vec;
};



class Matrix
{
public:

  // Constructors
  Matrix();
  Matrix(int, int);
  Matrix(G4double[]);
  
  // Destructor
  ~Matrix();

  inline int GetRows() const { return nr; }
  inline int GetCols() const { return nc; }

  // Puts control point into matrix location (i,j)
  inline void put(int i,int j, G4double x){ data[i*nc+j]=x; }
  
  // Retrieves control point from matrix location (i,j)
  inline G4double get(int i, int j) const
  {
    return data[i*nc+j];
  }


private:				      

  G4double* data;
  int nr, nc;
};


#endif







