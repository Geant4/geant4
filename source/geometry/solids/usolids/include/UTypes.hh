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
// UTypes
//
// Description:
//
//   Internal utility types defined for the unified solids library
//
// 19.10.12 Marek Gayer
// --------------------------------------------------------------------

#ifndef USOLIDS_Utypes
#define USOLIDS_Utypes

#include "UVector3.hh"
#include <iostream>
#include <string>
#include <vector>

class __void__;

typedef unsigned int UInt_t;

struct UBBoxStruct
{
  double extent[3];  // half-lengths on the 3 axis (arrays for indexing)
  double orig[3];    // center coordinates
};

/*struct UBuffer3DStruct {
   const int fType;          // Primitive type - predefined ones in TBuffer3DTypes.h

   UInt_t    fNbPnts;        // Number of points describing the shape
   UInt_t    fNbSegs;        // Number of segments describing the shape
   UInt_t    fNbPols;        // Number of polygons describing the shape

   UInt_t    fPntsCapacity;  // Current capacity of fPnts space
   UInt_t    fSegsCapacity;  // Current capacity of fSegs space
   UInt_t    fPolsCapacity;  // Current capacity of fSegs space

   UInt_t    fSections;      // Section validity flags
   double   *fPnts;              // x0, y0, z0, x1, y1, z1, ..... ..... ....
   int      *fSegs;              // c0, p0, q0, c1, p1, q1, ..... ..... ....
   int      *fPols;              // c0, n0, s0, s1, ... sn, c1, n1, s0, ... sn
};
*/
/*
struct UFacet2{
  UInt_t f1;//number of vertices from verticesList forming a facet
  UInt_t f2;
  UInt_t f3;
  UInt_t f4;
};
struct UPolyhedron2{
  std::vector<UVector3>  vertices;//List of Vertices for Polyhedron used in G4Vis
  std::vector<UFacet2>   facets;//List of Facets;
};
*/

typedef UBBoxStruct UBBox;
//typedef UBuffer3DStruct UBuffer3D;
typedef std::string UGeometryType;

#endif
