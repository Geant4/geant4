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
// $Id: G4ControlPoints.hh,v 1.6 2001-07-11 09:59:33 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ControlPoints
//
// Class Description:
//
// G4ControlPoints defines an array of G4PointRat which in turn
// are made of G4Point3D.

// Author: J.Sulkimo, P.Urban.
// Revisions by: A.Floquet, L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __G4ControlPoints_h
#define __G4ControlPoints_h 1

#include "G4PointRat.hh"

class G4ControlPoints 
{
public:  // with description

  G4ControlPoints();
  G4ControlPoints( G4int, G4int);
  G4ControlPoints( G4int, G4int, G4int );
  G4ControlPoints( const G4ControlPoints& );
    // Constructors.

  ~G4ControlPoints();
    // Destructor.

  G4ControlPoints& operator = ( const G4ControlPoints& );
    // Assignment operator.

  void SetWeights(G4double*);
   
  void CalcValues(G4double k1, G4double param, G4Point3D& pts1,
		  G4double k2, G4Point3D& pts2);
  void CalcValues(G4double k1, G4double param, G4PointRat& pts1,
		  G4double k2, G4PointRat& pts2);
  
  inline G4int GetRows() const;
  inline G4int GetCols() const;
  
  inline void put(G4int i, G4int j, const G4Point3D &tmp);
  inline void put(G4int i, G4int j, const G4PointRat& tmp);
    // Put control point into matrix location (i,j).

  inline G4Point3D Get3D(G4int i, G4int j) const;
  inline G4PointRat& GetRat(G4int i, G4int j) const;
    // Retrieve control point from matrix location (i,j);

  G4double ClosestDistanceToPoint(const G4Point3D&);


private:				      
    
  inline G4double Calc(G4double k1, G4double par, G4double old_val,
                       G4double k2, G4double new_val );

  G4PointRat** data;
    // G4PointRat datas are made of
    // - a point 3D
    // - an additional value: the scale factor which is set to 1 by default.

  G4int nr, nc;
};

#include "G4ControlPoints.icc"

#endif
