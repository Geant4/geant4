// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4KnotVector.hh,v 1.2 1999-11-08 09:50:31 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef __KNOTVECTOR_H
#define __KNOTVECTOR_H
#include <math.h>
//#include "STEPaggregate.h"
#include "geomdefs.hh"


class G4KnotVector
{
 public:
 
  G4KnotVector();
//  G4KnotVector(const int Size, const int* MultiList, STEPaggregate& Aggr);
//  G4KnotVector(const int Size, STEPaggregate& Aggr);
  G4KnotVector(const int sz);
  G4KnotVector(const G4KnotVector& old_kv);
  ~G4KnotVector();

  // Gets number of knots
  inline int GetSize()const {return k_size;};

  // Retrieves knot from knot vector index knot_number
  inline G4double GetKnot(const int knot_number){return knots[knot_number];}

  // Sets knot vector index knot_number to value
  inline void   PutKnot(const int knot_number, const G4double value)
  {
    knots[knot_number]=value;
  }

  // Adds the internal knots to the new knot vector
  G4KnotVector* MultiplyKnotVector( const int num, const G4double value);  

  // Creates the new vector by merging the old vector with the
  // knots in the vector knots_to_add
  G4double* MergeKnotVector( const G4double *knots_to_add, const int add_size);

  // Finds out how many Times val occurs in the knot vector
  int  CheckKnotVector(const G4double val);

  // Copies either the first half or the second half of
  // the new knot vector values to the knot vectors of the
  // new surfaces created by splitting
  void ExtractKnotVector( G4KnotVector* kv, const int upper, const int lower);

  // Searches the knot vector for the value and returns the index
  // This is used in the Evaluation of the intersection to find
  // out between which knots the intersection point is on the b-spline
  // surface
  int  GetKnotIndex(G4double k_value, const int order);


 private:

  // Number of knots
  int k_size;

  // Knot vector
  G4double *knots;

  inline G4double ApxEq(const G4double x,const G4double y) 
  {
    return (fabs(x - y) < kCarTolerance);
  }

};

#endif
