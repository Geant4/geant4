// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ControlPoints.hh,v 1.2 1999-12-15 14:49:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Modif 8 oct 98 : A.Floquet
//      G4PointRat datas are made of
// 	  . a point 3D
//	  . a additional value : the scale factor which is set to 1 by default
//
//      G4ControlPoints includes only G4PointRat which in turn are made 
//      of G4Point3D

#ifndef __G4ControlPoints_h
#define __G4ControlPoints_h 1

#include "G4PointRat.hh"

class G4ControlPoints 
{
public:

  // Constructors
  G4ControlPoints();

  G4ControlPoints(const STEPaggregate& Aggr, const int Rational);

  G4ControlPoints( int, int);

  G4ControlPoints( int , int , int );

  G4ControlPoints(const G4ControlPoints&);
  
  // Destructor
  ~G4ControlPoints();

  void SetWeights(G4double*);
   
  void CalcValues(G4double k1, G4double param, G4Point3D& pts1,
		  G4double k2, G4Point3D& pts2);
  void CalcValues(G4double k1, G4double param, G4PointRat& pts1,
		  G4double k2, G4PointRat& pts2);
  
  inline int GetRows() const {return nr;}
  inline int GetCols() const {return nc;}
  
  // Puts control point into matrix location (i,j)
  inline void put(const int i, const int j, const G4Point3D &tmp)
  {
    *data[i*nc+j]=tmp;	// tmp is converted to a PointRat
			// by the member affectation function
			// of the G4PointRat class
  }
  
  
  inline void put(const int i, const int j, const G4PointRat& tmp)
  {
    *data[i*nc+j]=tmp;
  }

  
  // Retrieves control point from matrix location (i,j)
  inline G4Point3D Get3D(const int i, const int j) const
  {
    return (data[i*nc+j])->pt();
  }
  
  inline G4PointRat& GetRat(const int i, const int j) const
  {
    return *data[i*nc+j];
  }        
      
  G4double ClosestDistanceToPoint(const G4Point3D&);


private:				      
    
  inline G4double Calc(const G4double k1, const G4double par, 
		       const G4double old_val, const G4double k2, 
		       const G4double new_val                     ) 
  {
    return (((k1 - par) * old_val +(par - k2) * new_val) / (k1-k2)); 
  }

  G4PointRat** data;
  int nr, nc;
};

#endif










