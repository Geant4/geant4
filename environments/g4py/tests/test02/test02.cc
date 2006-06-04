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
// $Id: test02.cc,v 1.3 2006-06-04 21:35:59 kmura Exp $
// ====================================================================
//   test02.cc
//
//   test of simple inheritance
//   - pure virtual class
//   - class inheritance
//   - function override
//   - polymorphic behavior
//   - pointer argument and down cast
//   - const flags
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

#include "XBase.hh"
#include "AClass.hh"
#include "BClass.hh"

using namespace boost::python;

BOOST_PYTHON_MODULE(test02) {
  class_<XBase, boost::noncopyable>("XBase", "Base Class", no_init)
    .add_property("ival", &XBase::GetIVal, &XBase::SetIVal)
    .add_property("dval", &XBase::GetDVal, &XBase::SetDVal)
    .def("AMethod", &XBase::AMethod)
    .def("VMethod", &XBase::VMethod)
    ;
  
  class_<AClass, bases<XBase> >( "AClass", "Derived Class A")
    .def(init<>())
    .def("AMethod", &AClass::AMethod)
    ;

  class_<BClass, bases<XBase> >( "BClass", "Derived Class B")
    .def(init<>())
    .def("AMethod", &BClass::AMethod)
    ;
}

