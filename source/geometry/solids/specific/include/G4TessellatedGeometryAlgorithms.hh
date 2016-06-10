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
// * technical work of the GEANT4 collaboration and of QinetiQ Ltd,   *
// * subject to DEFCON 705 IPR conditions.                            *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id:
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// Class G4TessellatedGeometryAlgorithms
//
// Class description:
//
//   The G4TessellatedGeometryAlgorithms class is used to contain standard
//   routines to determine whether (and if so where) simple geometric shapes
//   intersect.
//   
//   The constructor doesn't need to do anything, and neither does the
//   destructor.
//   
//   IntersectLineAndTriangle2D
//     Determines whether there is an intersection between a line defined
//     by r = p + s.v and a triangle defined by verticies P0, P0+E0 and P0+E1.
//     Here:
//        p = 2D vector
//        s = scaler on [0,infinity)
//        v = 2D vector
//        P0, E0 and E1 are 2D vectors
//     Information about where the intersection occurs is returned in the
//     variable location.
//
//   IntersectLineAndLineSegment2D
//     Determines whether there is an intersection between a line defined
//     by r = P0 + s.D0 and a line-segment with endpoints P1 and P1+D1.
//     Here:
//        P0 = 2D vector
//        s  = scaler on [0,infinity)
//        D0 = 2D vector
//        P1 and D1 are 2D vectors
//     Information about where the intersection occurs is returned in the
//     variable location.

// CHANGE HISTORY
// --------------
//
// 07 August 2007, P R Truscott, QinetiQ Ltd, UK - Created, with member
//                 functions based on the work of Rickard Holmberg.
// 12 October 2012, M Gayer, CERN, - Reviewed optimized implementation.
//
///////////////////////////////////////////////////////////////////////////////
#ifndef G4TessellatedGeometryAlgorithms_hh
#define G4TessellatedGeometryAlgorithms_hh 1

#include "G4TwoVector.hh"

class G4TessellatedGeometryAlgorithms
{
  public:

    static G4bool IntersectLineAndTriangle2D (const G4TwoVector &p,
                                              const G4TwoVector &v,
                                              const G4TwoVector &p0, 
                                              const G4TwoVector &e0,
                                              const G4TwoVector &e1,
                                                    G4TwoVector location[2]);
    static G4int IntersectLineAndLineSegment2D (const G4TwoVector &p0,
                                                const G4TwoVector &d0,
                                                const G4TwoVector &p1,
                                                const G4TwoVector &d1,
                                                      G4TwoVector location[2]);
    static G4double cross(const G4TwoVector &v1, const G4TwoVector &v2);
};

#endif
