// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ConvexHull.hh,v 2.1 1998/10/20 16:31:16 broglia Exp $
// GEANT4 tag $Name: geant4-00 $
//
#ifndef __CONVEXHULL_H
#define __CONVEXHULL_H


class G4ConvexHull
{
 public:

  G4ConvexHull *next;
  G4double      param;
  G4double      min;
  G4double      max;

  G4ConvexHull(){};
  G4ConvexHull(G4double pparam, G4double mmin, G4double mmax)
  { 
    next  = this;
    param = pparam;
    min   = mmin;
    max   = mmax;
  }
  
  ~G4ConvexHull(){}
};

#endif

