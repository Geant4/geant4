//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
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
