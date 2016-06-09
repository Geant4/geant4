// -*- C++ -*-
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of the output method of the HepRotation class,
// which was introduced when ZOOM PhysicsVectors was merged in.
//

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/Rotation.h"

#include <iomanip>
#include <iostream>

namespace CLHEP  {

std::ostream & HepRotation::print( std::ostream & os ) const {
  os << "\n   [ ( " <<
        std::setw(11) << std::setprecision(6) << xx() << "   " <<
        std::setw(11) << std::setprecision(6) << xy() << "   " <<
        std::setw(11) << std::setprecision(6) << xz() << ")\n"
     << "     ( " <<
        std::setw(11) << std::setprecision(6) << yx() << "   " <<
        std::setw(11) << std::setprecision(6) << yy() << "   " <<
        std::setw(11) << std::setprecision(6) << yz() << ")\n"
     << "     ( " <<
        std::setw(11) << std::setprecision(6) << zx() << "   " <<
        std::setw(11) << std::setprecision(6) << zy() << "   " <<
        std::setw(11) << std::setprecision(6) << zz() << ") ]\n";
	return os;
}


}  // namespace CLHEP
