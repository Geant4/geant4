// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ConvexHull.hh,v 1.3 2000-08-28 08:57:44 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ConvexHull
//
// Class Description:
//
// Definition of a convex hull list.

// Author: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __CONVEXHULL_H
#define __CONVEXHULL_H


class G4ConvexHull
{
 public:  // with description

  G4ConvexHull(){}
  G4ConvexHull(G4double pparam, G4double mmin, G4double mmax)
    : next(this), param (pparam), min(mmin), max(mmax){}
    // Constructors
  
  ~G4ConvexHull(){}
    // Destructor

 public:  // withouth description

  G4ConvexHull *next;
  G4double      param;
  G4double      min;
  G4double      max;
};

#endif
