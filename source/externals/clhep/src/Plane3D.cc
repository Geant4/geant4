// -*- C++ -*-
// $Id:$
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// Hep geometrical 3D Plane class
//
// Author: Evgeni Chernyaev <Evgueni.Tcherniaev@cern.ch>
//
// History:
// 22.09.96 E.Chernyaev - initial version
// 19.10.96 J.Allison - added == and <<.
// 15.04.03 E.Chernyaev - CLHEP-1.9: template version

#include <iostream>
#include "CLHEP/Geometry/Plane3D.h"

namespace HepGeom {
  //--------------------------------------------------------------------------
  std::ostream &
  operator<<(std::ostream & os, const Plane3D<float> & p) {
    return os
      << '(' << p.a() << ',' << p.b() << ',' << p.c() << ',' << p.d() << ')';
  }

  //--------------------------------------------------------------------------
  std::ostream &
  operator<<(std::ostream & os, const Plane3D<double> & p) {
    return os
      << '(' << p.a() << ',' << p.b() << ',' << p.c() << ',' << p.d() << ')';
  }
} /* namespace HepGeom */
