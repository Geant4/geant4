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
// $Id: test01.cc,v 1.3 2006-06-04 21:35:59 kmura Exp $
// ====================================================================
//   test01.cc
//
//   exposure of a simple class
//   - constructor
//   - multiple constructors
//   - default value in function arguments
//   - setter/getter
//   - method
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

#include "AClass.hh"

using namespace boost::python;

BOOST_PYTHON_MODULE(test01) {
  class_<AClass>( "AClass", "a simple class")
    .def(init<>())
    .def(init<int>())
    .def(init<int, double>())
    .add_property("ival", &AClass::GetIVal, &AClass::SetIVal)
    .add_property("dval", &AClass::GetDVal, &AClass::SetDVal)
    .def("AMethod", &AClass::AMethod)
    ;
}

