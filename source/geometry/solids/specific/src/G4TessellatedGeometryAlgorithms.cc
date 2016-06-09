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
//
// $Id: G4TessellatedGeometryAlgorithms.cc,v 1.6 2008-12-18 12:57:36 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4TessellatedGeometryAlgorithms.cc
//
// Date:                07/08/2005
// Author:              Rickard Holmberg & Pete Truscott
// Organisation:        QinetiQ Ltd, UK (PT)
// Customer:            ESA-ESTEC / TEC-EES
// Contract:            
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 07 August 2007, P R Truscott, QinetiQ Ltd, UK - Created, with member
//                 functions based on the work of Rickard Holmberg.
//
// 26 September 2007
//                 P R Truscott, qinetiQ Ltd, UK
//                 Updated to assign values of location array, not update
//                 just the pointer.
//
///////////////////////////////////////////////////////////////////////////////

#include "G4TessellatedGeometryAlgorithms.hh"
///////////////////////////////////////////////////////////////////////////////
//
// Pointer to single instance of class.
//
G4TessellatedGeometryAlgorithms* G4TessellatedGeometryAlgorithms::fInstance = 0;

///////////////////////////////////////////////////////////////////////////////
//
// G4TessellatedGeometryAlgorithms
//
// Constructor doesn't need to do anything since this class just allows access
// to the geometric algorithms contained in member functions.
//
G4TessellatedGeometryAlgorithms::G4TessellatedGeometryAlgorithms ()
{
}

///////////////////////////////////////////////////////////////////////////////
//
// GetInstance
//
// This is the access point for this singleton.
//
G4TessellatedGeometryAlgorithms* G4TessellatedGeometryAlgorithms::GetInstance()
{
  static G4TessellatedGeometryAlgorithms worldStdGeom;
  if (!fInstance)
  {
    fInstance = &worldStdGeom;
  }
  return fInstance;
}

///////////////////////////////////////////////////////////////////////////////
//
// IntersectLineAndTriangle2D
//
// Determines whether there is an intersection between a line defined
// by r = p + s.v and a triangle defined by verticies P0, P0+E0 and P0+E1.
//
// Here:
//        p = 2D vector
//        s = scaler on [0,infinity)
//        v = 2D vector
//        P0, E0 and E1 are 2D vectors
// Information about where the intersection occurs is returned in the
// variable location.
//
// This is based on the work of Rickard Holmberg.
//
G4bool G4TessellatedGeometryAlgorithms::IntersectLineAndTriangle2D (
  const G4TwoVector p,  const G4TwoVector v,
  const G4TwoVector P0, const G4TwoVector E0, const G4TwoVector E1,
  G4TwoVector location[2])
{
  G4TwoVector loc0[2];
  G4int e0i = IntersectLineAndLineSegment2D (p,v,P0,E0,loc0);
  if (e0i == 2)
  {
    location[0] = loc0[0];
    location[1] = loc0[1];
    return true;
  }
  
  G4TwoVector loc1[2];
  G4int e1i = IntersectLineAndLineSegment2D (p,v,P0,E1,loc1);
  if (e1i == 2)
  {
    location[0] = loc1[0];
    location[1] = loc1[1];
    return true;
  }
  
  if ((e0i == 1) && (e1i == 1))
  {
    if ((loc0[0]-p).mag2() < (loc1[0]-p).mag2())
    {
      location[0] = loc0[0];
      location[1] = loc1[0];
    }
    else
    {
      location[0] = loc1[0];
      location[1] = loc0[0];
    }
    return true;
  }
  
  G4TwoVector P1 = P0 + E0;
  G4TwoVector DE = E1 - E0;
  G4TwoVector loc2[2];
  G4int e2i = IntersectLineAndLineSegment2D (p,v,P1,DE,loc2);
  if (e2i == 2)
  {
    location[0] = loc2[0];
    location[1] = loc2[1];
    return true;
  }

  if ((e0i == 0) && (e1i == 0) && (e2i == 0)) { return false; }

  if ((e0i == 1) && (e2i == 1))
  {
    if ((loc0[0]-p).mag2() < (loc2[0]-p).mag2())
    {
      location[0] = loc0[0];
      location[1] = loc2[0];
    }
    else
    {
      location[0] = loc2[0];
      location[1] = loc0[0];
    }
    return true;
  }

  if ((e1i == 1) && (e2i == 1))
  {
    if ((loc1[0]-p).mag2() < (loc2[0]-p).mag2())
    {
      location[0] = loc1[0];
      location[1] = loc2[0];
    }
    else
    {
      location[0] = loc2[0];
      location[1] = loc1[0];
    }
    return true;
  }

  return false;
}

///////////////////////////////////////////////////////////////////////////////
//
// IntersectLineAndLineSegment2D
//
// Determines whether there is an intersection between a line defined
// by r = P0 + s.D0 and a line-segment with endpoints P1 and P1+D1.
// Here:
//        P0 = 2D vector
//        s  = scaler on [0,infinity)
//        D0 = 2D vector
//        P1 and D1 are 2D vectors
//
// This function returns:
// 0 - if there is no intersection;
// 1 - if there is a unique intersection;
// 2 - if the line and line-segments overlap, and the intersection is a
//     segment itself.
// Information about where the intersection occurs is returned in the
// as ??.
//
// This is based on the work of Rickard Holmberg as well as material published
// by Philip J Schneider and David H Eberly, "Geometric Tools for Computer
// Graphics," ISBN 1-55860-694-0, pp 244-245, 2003.
//
G4int G4TessellatedGeometryAlgorithms::IntersectLineAndLineSegment2D (
  const G4TwoVector P0, const G4TwoVector D0,
  const G4TwoVector P1, const G4TwoVector D1,
  G4TwoVector location[2])
{
  G4TwoVector E     = P1 - P0;
  G4double kross    = cross(D0,D1);
  G4double sqrKross = kross * kross;
  G4double sqrLen0  = D0.mag2();
  G4double sqrLen1  = D1.mag2();
  location[0]       = G4TwoVector(0.0,0.0);
  location[1]       = G4TwoVector(0.0,0.0);

  if (sqrKross > DBL_EPSILON * DBL_EPSILON * sqrLen0 * sqrLen1)
  {
//
//
// The line and line segment are not parallel.  Determine if the intersection
// is in positive s where r = P0 + s*D0, and for 0<=t<=1 where r = p1 + t*D1.
//
    G4double s = cross(E,D1)/kross;
    if (s < 0)          return 0; // Intersection does not occur for positive s.
    G4double t = cross(E,D0)/kross;
    if (t < 0 || t > 1) return 0; // Intersection does not occur on line-segment.
//
//
// Intersection of lines is a single point on the forward-propagating line
// defined by r = P0 + s*D0, and the line segment defined by  r = p1 + t*D1.
//
    location[0] = P0 + s*D0;
    return 1;
  }
//
//
// Line and line segment are parallel.  Determine whether they overlap or not.
//
  G4double sqrLenE = E.mag2();
  kross            = cross(E,D0);
  sqrKross         = kross * kross;
  if (sqrKross > DBL_EPSILON * DBL_EPSILON * sqrLen0 * sqrLenE)
  {
    return 0; //Lines are different.
  }
//
//
// Lines are the same.  Test for overlap.
//
  G4double s0   = D0.dot(E)/sqrLen0;
  G4double s1   = s0 + D0.dot(D1)/sqrLen0;
  G4double smin = 0.0;
  G4double smax = 0.0;

  if (s0 < s1) {smin = s0; smax = s1;}
  else         {smin = s1; smax = s0;}

  if (smax < 0.0) return 0;
  else if (smin < 0.0)
  {
    location[0] = P0;
    location[1] = P0 + smax*D0;
    return 2;
  }
  else
  {
    location[0] = P0 + smin*D0;
    location[1] = P0 + smax*D0;
    return 2;
  }
}
