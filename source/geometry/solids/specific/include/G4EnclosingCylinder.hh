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
// $Id: G4EnclosingCylinder.hh,v 1.5 2002-10-28 11:47:50 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4EnclosingCylinder
//
// Class description:
//
//   Definition of a utility class for quickly deciding if a point
//   is clearly outside a polyhedra or polycone or deciding if
//   a trajectory is clearly going to miss those shapes.

// Author: 
//   David C. Williams (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#ifndef G4EnclosingCylinder_hh
#define G4EnclosingCylinder_hh

#include "globals.hh"
#include "geomdefs.hh"
#include "G4ThreeVector.hh"

class G4ReduciblePolygon;

class G4EnclosingCylinder
{
  public:  // with description

  G4EnclosingCylinder( const G4ReduciblePolygon *rz,
                             G4bool phiIsOpen, 
                             G4double startPhi, G4double totalPhi );
  ~G4EnclosingCylinder();
  
  G4bool MustBeOutside( const G4ThreeVector &p ) const;
    // Decide very rapidly if the point is outside the cylinder.
    // If one is not certain, return false.

  G4bool ShouldMiss( const G4ThreeVector &p, const G4ThreeVector &v ) const;
    // Decide very rapidly if the trajectory is going to miss the cylinder.
    // If one is not sure, return false.

  protected:

  G4double radius;    // radius of our cylinder
  G4double zLo, zHi;  // z extent
  
  G4bool    phiIsOpen; // true if there is a phi segment
  G4double  startPhi,  // for isPhiOpen==true, starting of phi segment
            totalPhi;  // for isPhiOpen==true, size of phi segment

  G4double rx1, ry1,
           dx1, dy1;
  G4double rx2, ry2,
           dx2, dy2;
     
  G4bool   concave;  // true, if x/y cross section is concave
    
};

#endif
