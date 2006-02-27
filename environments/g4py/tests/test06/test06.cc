// $Id: test06.cc,v 1.1 2006-02-27 10:05:25 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   test06.cc
//
//   test for inheritance
//   - inherit C++ base class in Python side
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

#include "XBase.hh"
#include "ZBase.hh"

using namespace boost::python;

struct w_XBase : XBase, wrapper<XBase> {
  std::string PVMethod() {
    return this-> get_override("PVMethod")();
  }
};


struct w_ZBase : ZBase, wrapper<ZBase> {
  void VMethod(std::string message) {
    if(override f= this-> get_override("VMethod"))
      f(message);
    else 
      ZBase::VMethod(message);
  }

  void d_VMethod(std::string message) {
    this-> ZBase::VMethod(message);
  }
};


BOOST_PYTHON_MODULE(test06) {
  //class_<w_XBase, boost::noncopyable>("XBase", "base class", no_init)
  class_<w_XBase, boost::noncopyable>("XBase", "base class")
    .add_property("ival", &XBase::GetIVal, &XBase::SetIVal)
    .def("PVMethod", pure_virtual(&XBase::PVMethod))
    ;

  class_<w_ZBase, boost::noncopyable>("ZBase", "base class")
    .def("AMethod", &ZBase::AMethod)
    .def("VMethod", &ZBase::VMethod, &w_ZBase::d_VMethod)
    ;
}
