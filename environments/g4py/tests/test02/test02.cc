// $Id: test02.cc,v 1.1 2006-02-27 10:08:30 kmura Exp $
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

