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
// $Id: G4ConvexHull.hh,v 1.6 2002-02-26 11:16:55 gcosmo Exp $
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
    : param (pparam), min(mmin), max(mmax) { next=this; }
    // Constructors
  
  ~G4ConvexHull(){}
    // Destructor

  G4ConvexHull* GetNextHull() { return next; }
  G4double GetParam() const { return param; }
  G4double GetMin() const { return min; }
  G4double GetMax() const { return max; }

  void SetNextHull(G4ConvexHull* n) { next=n; }
  void SetParam(G4double p) { param=p; }
  void SetMin(G4double x) { min=x; }
  void SetMax(G4double y) { max=y; }

 private:

  G4ConvexHull *next;
  G4double      param;
  G4double      min;
  G4double      max;
};

#endif
