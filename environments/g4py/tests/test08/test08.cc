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
// $Id: test08.cc,v 1.3 2006-06-04 21:36:00 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   test08.cc
//
//   test for static member function
//
//                                         2005 Q
// ====================================================================

class AClass {
public:
  AClass() { }
  ~AClass() { }
  static int AMethod() { return 1; } 

};

// Boost.Python...
#include <boost/python.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(test08)
{
  class_<AClass>( "AClass", "a class")
    .def(init<>())
    .def("AMethod", &AClass::AMethod)
    .staticmethod("AMethod")
    ;
}

