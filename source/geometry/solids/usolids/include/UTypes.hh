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

typedef UBBoxStruct UBBox;
typedef std::string UGeometryType;

#endif
