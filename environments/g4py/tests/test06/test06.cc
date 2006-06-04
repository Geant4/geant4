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
// $Id: test06.cc,v 1.3 2006-06-04 21:35:59 kmura Exp $
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
