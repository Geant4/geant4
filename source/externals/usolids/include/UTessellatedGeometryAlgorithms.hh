//
// ********************************************************************
// * This Software is part of the AIDA Unified Solids Library package *
// * See: https://aidasoft.web.cern.ch/USolids                        *
// ********************************************************************
//
// $Id:$
//
// --------------------------------------------------------------------
//
// UTessellatedGeometryAlgorithms
//
// Class description:
//
// The UTessellatedGeometryAlgorithms class is used to contain standard
// routines to determine whether (and if so where) simple geometric shapes
// intersect.
//
// The constructor doesn't need to do anything, and neither does the
// destructor.
//
// IntersectLineAndTriangle2D
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
// IntersectLineAndLineSegment2D
//     Determines whether there is an intersection between a line defined
//     by r = P0 + s.D0 and a line-segment with endpoints P1 and P1+D1.
//     Here:
//        P0 = 2D vector
//        s  = scaler on [0,infinity)
//        D0 = 2D vector
//        P1 and D1 are 2D vectors
//     Information about where the intersection occurs is returned in the
//     variable location.
//
// 11.07.12 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#ifndef UTessellatedGeometryAlgorithms_hh
#define UTessellatedGeometryAlgorithms_hh 1

#include "UVector2.hh"

class UTessellatedGeometryAlgorithms
{
  public:

    static bool IntersectLineAndTriangle2D(const UVector2& p,
                                           const UVector2& v,
                                           const UVector2& p0,
                                           const UVector2& e0,
                                           const UVector2& e1,
                                           UVector2 location[2]);

    static int IntersectLineAndLineSegment2D(const UVector2& p0,
                                             const UVector2& d0,
                                             const UVector2& p1,
                                             const UVector2& d1,
                                             UVector2 location[2]);

    static double Cross(const UVector2& v1, const UVector2& v2);

};

#endif
